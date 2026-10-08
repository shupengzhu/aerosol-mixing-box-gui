"""Local web interface: python -m scrambox serve [--port 8765] [--host 127.0.0.1]

Standard-library HTTP server with a JSON API and a static single-page front end
(scrambox/web). On a shared server every user starts their own instance (their runs and
saved scenarios live in SCRAMBOX_USER_DIR, default ~/.scrambox) and reaches it through an
ssh tunnel:  ssh -L 8765:localhost:8765 <server>.

API
  GET    /api/meta                     enums, defaults, species tables, kernels
  GET    /api/scenarios                scenario index;  GET /api/scenarios/<id>  full scenario
  POST   /api/scenarios                save a user scenario;  DELETE /api/scenarios/<id>
  POST   /api/describe                 information panel for a case
  POST   /api/grid                     sections, classes and labels for a case
  POST   /api/runs                     start a run {case, name?} -> {id}
  GET    /api/runs                     list;  GET /api/runs/<id>  status;  DELETE /api/runs/<id>
  GET    /api/runs/<id>/results        arrays for the result views (+ post-processing)
  POST   /api/runs/<id>/post           (re)run post-processing with options
  GET    /api/runs/<id>/case           the resolved case of a run
  POST   /api/experiments              start an experiment set;  GET /api/experiments
  GET    /api/experiments/<id>         status and comparison metrics
"""
from __future__ import annotations

import json
import mimetypes
import os
import re
import shutil
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np

from . import case as casemod
from . import scenarios as scn
from . import species as spmod

WEB = Path(__file__).resolve().parent / "web"
POOL = ThreadPoolExecutor(max_workers=int(os.environ.get("SCRAMBOX_WORKERS", "4")))
LOCK = threading.Lock()


def user_dir() -> Path:
    d = Path(os.environ.get("SCRAMBOX_USER_DIR", Path.home() / ".scrambox"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def runs_dir() -> Path:
    d = user_dir() / "runs"
    d.mkdir(exist_ok=True)
    return d


def exps_dir() -> Path:
    d = user_dir() / "experiments"
    d.mkdir(exist_ok=True)
    return d


def _safe_id(s: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", s or ""):
        raise KeyError(s)
    return s


def _new_id(name: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", name or "run").strip("_")[:32] or "run"
    return time.strftime("%Y%m%d_%H%M%S_") + slug + "_" + uuid.uuid4().hex[:4]


def _round(a, sig=5):
    a = np.asarray(a, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag = np.where(a != 0, np.floor(np.log10(np.abs(a))), 0)
    f = 10.0 ** (sig - 1 - mag)
    return (np.round(a * f) / f).tolist()


# ------------------------------------------------------------------ meta
def meta() -> dict:
    from .experiments import PRESET_DIMENSIONS
    from .grid import PRESETS
    from .post import ccn, hetchem, optics
    from .post.common import KAPPA
    from .run import KERNEL_BUILD
    libs = {}
    for prof in ("iso", "noiso"):
        t = spmod.library(prof)
        d = t.to_dict()
        d["kappa"] = {s: KAPPA.get(s) for s in t.aerosol}
        libs[prof] = d
    kernels = {}
    for v in ("v1", "v23"):
        f = KERNEL_BUILD / v / "KERNEL_VERSION"
        kernels[v] = {"available": (KERNEL_BUILD / v / "scram_box").exists(),
                      "version": f.read_text() if f.exists() else ""}
    return {"defaults": casemod.DEFAULTS, "enums": casemod.ENUMS, "species": libs,
            "kernels": kernels, "categories": scn.CATEGORIES,
            "mixing_presets": {k: {"external": v[0], "primary": v[1], "secondary": v[2]}
                               for k, v in PRESETS.items()},
            "post_defaults": {"optics": optics.DEFAULTS, "ccn": ccn.DEFAULTS,
                              "hetchem": hetchem.DEFAULTS, "gamma_cesm": hetchem.GAMMA_CESM},
            "experiment_dimensions": PRESET_DIMENSIONS}


def grid_info(case: dict) -> dict:
    from .compile import build_grid
    from .grid import class_label
    c = casemod.resolve(case)
    g = build_grid(c)
    t = g["table"]
    return {"bounds_um": g["bounds_um"], "groups": t.groups,
            "classes": [[list(b) for b in cls] for cls in g["classes"]],
            "class_labels": [class_label(cls, t.groups) for cls in g["classes"]],
            "aerosol": t.aerosol, "gases": t.gases, "group_of": t.group_of,
            "substep_s": c["time"]["substep_s"]}


# ------------------------------------------------------------------ runs
def _status(run: Path) -> dict:
    f = run / "status.json"
    st = json.loads(f.read_text()) if f.exists() else {"state": "unknown"}
    if st.get("state") == "running":
        out = run / "output.txt"
        try:
            n_done = sum(1 for line in open(out) if line.startswith("time")) if out.exists() else 0
            c = json.loads((run / "case.resolved.json").read_text())
            n_tot = int(c["time"]["duration_s"] // c["time"]["output_every_s"]) + 1
            st["progress"] = min(n_done / max(n_tot, 1), 1.0)
        except Exception:                       # noqa: BLE001
            pass
    return st


def _write_status(run: Path, **kw):
    f = run / "status.json"
    st = json.loads(f.read_text()) if f.exists() else {}
    st.update(kw)
    f.write_text(json.dumps(st, indent=1))


def start_run(case: dict, name: str | None = None) -> str:
    casemod.resolve(case)                        # validate before queueing
    rid = _new_id(name or case.get("case", {}).get("name", "run"))
    run = runs_dir() / rid
    run.mkdir()
    (run / "case.json").write_text(json.dumps(case, indent=1, ensure_ascii=False))
    _write_status(run, state="queued", name=name or case.get("case", {}).get("name", ""), created=time.time())

    def job():
        from .run import run_case
        _write_status(run, state="running", started=time.time())
        try:
            out = run_case(case, run)
            _write_status(run, state="done", finished=time.time(), cpu_s=out["summary"]["cpu_s"])
        except Exception as e:                  # noqa: BLE001
            _write_status(run, state="error", finished=time.time(), error=str(e)[-1500:])
    POOL.submit(job)
    return rid


def list_runs() -> list[dict]:
    out = []
    for d in sorted(runs_dir().iterdir(), reverse=True):
        if d.is_dir():
            st = _status(d)
            out.append({"id": d.name, **{k: st.get(k) for k in ("state", "name", "created", "cpu_s", "progress", "error")}})
    return out


def run_results(rid: str) -> dict:
    from .run import load_run
    run = runs_dir() / _safe_id(rid)
    r = load_run(run)
    res, g = r["results"], r["grid"]
    return {"id": rid, "name": r["case"]["case"]["name"], "summary": r["summary"],
            "time_s": res["time_s"].tolist(), "met": _round(res["met"]),
            "bounds_um": g["bounds_um"], "class_labels": g["class_labels"],
            "species": g["species"]["aerosol"], "gases": g["species"]["gases"],
            "number_m3": _round(res["number_m3"]), "mass_ugm3": _round(res["mass_ugm3"]),
            "gas_ugm3": _round(res["gas_ugm3"]), "post": r.get("post"),
            "het": _round(res["het"]) if "het" in res else None,
            "kernel": r["case"]["kernel"]["isorropia"], "notes": g.get("notes", [])}


def run_post(rid: str, cfg: dict) -> dict:
    from .post import run_all
    from .run import load_run
    run = runs_dir() / _safe_id(rid)
    r = load_run(run)
    return run_all(r, cfg, run)


# ------------------------------------------------------------------ experiments
def start_experiment(spec: dict) -> str:
    from .experiments import expand, run_set
    for c in expand(spec):
        casemod.resolve(c)
    eid = _new_id(spec.get("name", "experiment"))
    d = exps_dir() / eid
    d.mkdir()
    (d / "status.json").write_text(json.dumps({"state": "queued", "variants": []}))
    POOL.submit(lambda: run_set(spec, d, workers=int(os.environ.get("SCRAMBOX_EXP_WORKERS", "4"))))
    return eid


def list_experiments() -> list[dict]:
    out = []
    for d in sorted(exps_dir().iterdir(), reverse=True):
        if not d.is_dir():
            continue
        st = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
        spec = json.loads((d / "experiment.json").read_text()) if (d / "experiment.json").exists() else {}
        out.append({"id": d.name, "name": spec.get("name"), "state": st.get("state"),
                    "variants": [v.get("label") for v in spec.get("variants", [])]})
    return out


def experiment(eid: str) -> dict:
    from .experiments import compare
    d = exps_dir() / _safe_id(eid)
    st = json.loads((d / "status.json").read_text())
    if st.get("state") in ("queued", "running") or not (d / "experiment.json").exists():
        return {"id": eid, "status": st.get("state"), "progress": st.get("variants", [])}
    out = compare(d)
    out["id"] = eid
    return out


# ------------------------------------------------------------------ HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = "scrambox/0.2"

    def log_message(self, fmt, *args):         # quiet
        pass

    def _send(self, code, body, ctype="application/json; charset=utf-8"):
        data = body if isinstance(body, bytes) else json.dumps(body, ensure_ascii=False, allow_nan=False,
                                                             default=_json_default).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        self.wfile.write(data)

    def _body(self):
        n = int(self.headers.get("Content-Length") or 0)
        return json.loads(self.rfile.read(n) or b"{}")

    def _route(self, method):
        path = unquote(urlparse(self.path).path)
        parts = [p for p in path.split("/") if p]
        try:
            if not parts or parts[0] != "api":
                return self._static(path)
            p = parts[1:]
            if method == "GET":
                if p == ["meta"]:
                    return self._send(200, meta())
                if p == ["scenarios"]:
                    return self._send(200, scn.index())
                if len(p) == 2 and p[0] == "scenarios":
                    return self._send(200, scn.get(p[1]))
                if p == ["runs"]:
                    return self._send(200, list_runs())
                if len(p) == 2 and p[0] == "runs":
                    return self._send(200, _status(runs_dir() / _safe_id(p[1])))
                if len(p) == 3 and p[0] == "runs" and p[2] == "results":
                    return self._send(200, run_results(p[1]))
                if len(p) == 3 and p[0] == "runs" and p[2] == "case":
                    d = runs_dir() / _safe_id(p[1])
                    f = d / "case.json"
                    return self._send(200, json.loads(f.read_text()))
                if p == ["experiments"]:
                    return self._send(200, list_experiments())
                if len(p) == 2 and p[0] == "experiments":
                    return self._send(200, experiment(p[1]))
            elif method == "POST":
                body = self._body()
                if p == ["describe"]:
                    return self._send(200, scn.describe(body))
                if p == ["grid"]:
                    return self._send(200, grid_info(body))
                if p == ["resolve"]:
                    return self._send(200, casemod.resolve(body))
                if p == ["scenarios"]:
                    return self._send(200, scn.save_user(body))
                if p == ["runs"]:
                    return self._send(200, {"id": start_run(body["case"], body.get("name"))})
                if len(p) == 3 and p[0] == "runs" and p[2] == "post":
                    return self._send(200, run_post(p[1], body))
                if p == ["experiments"]:
                    return self._send(200, {"id": start_experiment(body)})
            elif method == "DELETE":
                if len(p) == 2 and p[0] == "scenarios":
                    scn.delete_user(p[1])
                    return self._send(200, {"deleted": p[1]})
                if len(p) == 2 and p[0] in ("runs", "experiments"):
                    base = runs_dir() if p[0] == "runs" else exps_dir()
                    d = base / _safe_id(p[1])
                    if d.exists():
                        shutil.rmtree(d)
                    return self._send(200, {"deleted": p[1]})
            return self._send(404, {"error": "not found"})
        except KeyError as e:
            return self._send(404, {"error": f"not found: {e}"})
        except (ValueError, FileNotFoundError) as e:
            return self._send(400, {"error": str(e)})
        except Exception as e:                  # noqa: BLE001
            traceback.print_exc()
            return self._send(500, {"error": f"{type(e).__name__}: {e}"})

    def _static(self, path):
        rel = "index.html" if path in ("", "/") else path.lstrip("/")
        f = (WEB / rel).resolve()
        if WEB not in f.parents or not f.is_file():
            return self._send(404, {"error": "not found"})
        ctype = mimetypes.guess_type(f.name)[0] or "application/octet-stream"
        if ctype.startswith("text/") or ctype in ("application/javascript",):
            ctype += "; charset=utf-8"
        return self._send(200, f.read_bytes(), ctype)

    def do_GET(self):
        self._route("GET")

    def do_POST(self):
        self._route("POST")

    def do_DELETE(self):
        self._route("DELETE")


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    raise TypeError(type(o).__name__)


def serve(host: str = "127.0.0.1", port: int = 8765):
    httpd = ThreadingHTTPServer((host, port), Handler)
    print(f"SCRAM box app on http://{host}:{port}  (user data: {user_dir()})")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
