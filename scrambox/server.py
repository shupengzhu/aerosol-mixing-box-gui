"""Web interface server: python -m scrambox serve [--port 8765] [--host 127.0.0.1]

Standard-library HTTP server with a JSON API and a static single-page front end
(scrambox/web). Two modes:

  local   no access code configured: one user, data in SCRAMBOX_USER_DIR (~/.scrambox)
  shared  an access code is set (SCRAMBOX_ACCESS_CODE or the file <data>/access_code):
          users sign in with the code and a name; each name gets its own workspace
          <data>/users/<id>/ (runs, experiment sets, saved scenarios); compute limits apply.
          Meant to sit behind an HTTPS reverse proxy (e.g. tailscale funnel) on 127.0.0.1.

API (all under /api, JSON)
  GET  session | meta | scenarios | scenarios/<id> | runs | runs/<id> | runs/<id>/results
       runs/<id>/case | runs/<id>/export | experiments | experiments/<id>
  POST login {code, name} | logout | describe | grid | resolve | scenarios | runs {case, name}
       runs/<id>/post | experiments
  DELETE scenarios/<id> | runs/<id> | experiments/<id>
"""
from __future__ import annotations

import hashlib
import hmac
import json
import mimetypes
import os
import re
import secrets
import shutil
import stat
import threading
import time
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from http.cookies import SimpleCookie
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import unquote, urlparse

import numpy as np

from . import case as casemod
from . import scenarios as scn
from . import species as spmod

WEB = Path(__file__).resolve().parent / "web"
POOL = ThreadPoolExecutor(max_workers=int(os.environ.get("SCRAMBOX_WORKERS", "6")))
LOCK = threading.Lock()

# limits (shared mode)
MAX_BODY = 8 * 1024 * 1024
MAX_EXPORT_BODY = 128 * 1024 * 1024
MAX_EXPORT_FILES = 1024
MAX_ACTIVE_PER_USER = int(os.environ.get("SCRAMBOX_MAX_ACTIVE", "3"))
MAX_COST = float(os.environ.get("SCRAMBOX_MAX_COST", "2e7"))      # kernel calls x sections
MAX_VARIANTS = int(os.environ.get("SCRAMBOX_MAX_VARIANTS", "10"))
RUN_TIMEOUT_S = float(os.environ.get("SCRAMBOX_RUN_TIMEOUT", "3600"))
MAX_PARTMC_PART = int(os.environ.get("SCRAMBOX_MAX_PARTMC_PART", "20000"))
MAX_PARTMC_COST = float(os.environ.get("SCRAMBOX_MAX_PARTMC_COST", "3e7"))   # particles x steps
PARTMC_TIMEOUT_S = float(os.environ.get("SCRAMBOX_PARTMC_TIMEOUT", "7200"))
SESSION_DAYS = 30


# ------------------------------------------------------------------ data locations
def data_dir() -> Path:
    d = Path(os.environ.get("SCRAMBOX_USER_DIR", Path.home() / ".scrambox"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def access_code() -> str | None:
    code = os.environ.get("SCRAMBOX_ACCESS_CODE")
    f = data_dir() / "access_code"
    if not code and f.exists():
        code = f.read_text().strip()
    return code or None


def shared_mode() -> bool:
    return access_code() is not None


class Workspace:
    """Directory tree of one user (or of the single local user)."""

    def __init__(self, root: Path, user: str | None = None):
        self.root, self.user = root, user
        root.mkdir(parents=True, exist_ok=True)

    @property
    def runs(self) -> Path:
        d = self.root / "runs"
        d.mkdir(exist_ok=True)
        return d

    @property
    def exps(self) -> Path:
        d = self.root / "experiments"
        d.mkdir(exist_ok=True)
        return d


def user_slug(name: str) -> str:
    ascii_part = re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")[:24]
    return (ascii_part + "_" if ascii_part else "u_") + hashlib.sha1(name.encode("utf-8")).hexdigest()[:8]


# ------------------------------------------------------------------ sessions
_SESS_LOCK = threading.Lock()
_FAILS: dict[str, list[float]] = {}


def _sessions_file() -> Path:
    return data_dir() / "sessions.json"


def _load_sessions() -> dict:
    f = _sessions_file()
    try:
        return json.loads(f.read_text()) if f.exists() else {}
    except Exception:                                   # noqa: BLE001
        return {}


def _save_sessions(s: dict) -> None:
    f = _sessions_file()
    tmp = f.with_suffix(".tmp")
    tmp.write_text(json.dumps(s))
    os.chmod(tmp, 0o600)
    tmp.replace(f)


def login(code: str, name: str, ip: str) -> str:
    now = time.time()
    with _SESS_LOCK:
        fails = [t for t in _FAILS.get(ip, []) if now - t < 600]
        if len(fails) >= 10:
            raise PermissionError("too many failed attempts, try again in 10 minutes")
        ok = hmac.compare_digest((code or "").strip().encode(), (access_code() or "").encode())
        name = (name or "").strip()
        if not ok:
            _FAILS[ip] = fails + [now]
            raise PermissionError("wrong access code")
        if not (1 <= len(name) <= 32) or not re.fullmatch(r"[\w .@-]+", name):
            raise ValueError("name: 1-32 letters, digits, Chinese characters, space . @ - _")
        token = secrets.token_hex(24)
        s = {k: v for k, v in _load_sessions().items() if now - v.get("t", 0) < SESSION_DAYS * 86400}
        s[token] = {"user": name, "t": now}
        _save_sessions(s)
    ws = workspace_for(name)
    prof = ws.root / "profile.json"
    if not prof.exists():
        prof.write_text(json.dumps({"name": name, "created": now}, ensure_ascii=False))
    return token


def logout(token: str | None) -> None:
    with _SESS_LOCK:
        s = _load_sessions()
        if token in s:
            del s[token]
            _save_sessions(s)


def session_user(token: str | None) -> str | None:
    if not token:
        return None
    v = _load_sessions().get(token)
    if not v or time.time() - v.get("t", 0) > SESSION_DAYS * 86400:
        return None
    return v["user"]


def workspace_for(user: str | None) -> Workspace:
    if user is None:
        return Workspace(data_dir())
    if re.fullmatch(r"5002:[1-9][0-9]*", user):
        # Signed website identities use a separate private tenant namespace.
        site_root = data_dir() / "site_users" / "5002" / user.split(":", 1)[1]
        site_root.mkdir(mode=0o700, parents=True, exist_ok=True)
        site_root.chmod(0o700)
        return Workspace(site_root, user)
    return Workspace(data_dir() / "users" / user_slug(user), user)


# ------------------------------------------------------------------ helpers
def _safe_id(s: str) -> str:
    if not re.fullmatch(r"[A-Za-z0-9_.-]{1,80}", s or "") or s.startswith("."):
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


def check_cost(case: dict) -> dict:
    """Validate a case and estimate its kernel cost; raise ValueError beyond the limit."""
    from .compile import build_grid
    c = casemod.resolve(case)
    g = build_grid(c)
    sections = (len(g["bounds_um"]) - 1) * len(g["classes"])
    calls = c["time"]["duration_s"] / c["time"]["substep_s"]
    cost = calls * sections
    if shared_mode() and cost > MAX_COST:
        raise ValueError(f"run too large for the shared server: {calls:.0f} kernel calls x {sections} sections "
                         f"= {cost:.2e} > {MAX_COST:.0e}; shorten the run, enlarge the sub-step or use fewer "
                         f"bins/classes (or install the standalone package)")
    return {"calls": calls, "sections": sections, "cost": cost}


def partmc_options(case: dict, opts: dict | None) -> dict:
    """Validated PartMC options for a case (shared-server limits apply)."""
    from .partmc import DEFAULTS, available
    o = dict(DEFAULTS, **(opts or {}))
    o = {"n_part": int(o["n_part"]), "del_t": float(o["del_t"]), "seed": int(o["seed"]),
         "weighting_exponent": float(o["weighting_exponent"])}
    c = casemod.resolve(case)
    if c["kernel"]["isorropia"] == "none":
        raise ValueError("PartMC runs use SCRAM's ISORROPIA scheme: choose ISORROPIA v2.3 or v1")
    if not available().get("v23" if c["kernel"]["isorropia"] == "v23" else "v1"):
        raise ValueError("PartMC is not installed on this server")
    if not (100 <= o["n_part"] <= 200000) or not (1.0 <= o["del_t"] <= 600.0):
        raise ValueError("PartMC: n_part 100-200000, del_t 1-600 s")
    if shared_mode():
        steps = c["time"]["duration_s"] / o["del_t"]
        if o["n_part"] > MAX_PARTMC_PART or o["n_part"] * steps > MAX_PARTMC_COST:
            raise ValueError(f"PartMC run too large for the shared server (at most {MAX_PARTMC_PART} particles and "
                             f"particles x steps <= {MAX_PARTMC_COST:.0e}; now {o['n_part']} x {steps:.0f})")
    return o


def _active_jobs(ws: Workspace) -> int:
    n = 0
    for base in (ws.runs, ws.exps):
        for d in base.iterdir():
            f = d / "status.json"
            if f.exists():
                try:
                    if json.loads(f.read_text()).get("state") in ("queued", "running"):
                        n += 1
                except Exception:                       # noqa: BLE001
                    pass
    return n


def _check_quota(ws: Workspace):
    if shared_mode() and _active_jobs(ws) >= MAX_ACTIVE_PER_USER:
        raise ValueError(f"at most {MAX_ACTIVE_PER_USER} runs or experiment sets at a time; wait for one to finish")


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
            "experiment_dimensions": PRESET_DIMENSIONS,
            "partmc": {"available": _partmc_available(), "defaults": _partmc_defaults()},
            "limits": ({"max_cost": MAX_COST, "max_active": MAX_ACTIVE_PER_USER, "max_variants": MAX_VARIANTS,
                        "max_partmc_part": MAX_PARTMC_PART, "max_partmc_cost": MAX_PARTMC_COST}
                       if shared_mode() else None)}


def _partmc_available():
    from .partmc import available
    return available()


def _partmc_defaults():
    from .partmc import DEFAULTS
    return DEFAULTS


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
    if not run.is_dir():
        raise KeyError(run.name)
    st = json.loads(f.read_text()) if f.exists() else {"state": "unknown"}
    if st.get("state") == "running" and st.get("model") == "partmc":
        from .partmc import progress
        st["progress"] = progress(run)
    elif st.get("state") == "running":
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


def start_run(ws: Workspace, case: dict, name: str | None = None, model: str = "scram",
              pm_opts: dict | None = None) -> str:
    if model == "partmc":
        pm = partmc_options(case, pm_opts)
    else:
        check_cost(case)
    _check_quota(ws)
    rid = _new_id(("PartMC " if model == "partmc" else "") + (name or case.get("case", {}).get("name", "run")))
    run = ws.runs / rid
    run.mkdir()
    (run / "case.json").write_text(json.dumps(case, indent=1, ensure_ascii=False))
    _write_status(run, state="queued", name=name or case.get("case", {}).get("name", ""), created=time.time(),
                  model=model)

    def job():
        _write_status(run, state="running", started=time.time())
        try:
            if model == "partmc":
                from .partmc import run_partmc
                t0 = time.time()
                run_partmc(case, run, n_part=pm["n_part"], del_t=pm["del_t"], seed=pm["seed"],
                           weighting_exponent=pm["weighting_exponent"], timeout_s=PARTMC_TIMEOUT_S, post=True)
                _write_status(run, state="done", finished=time.time(), cpu_s=time.time() - t0)
            else:
                from .run import run_case
                out = run_case(case, run, timeout_s=RUN_TIMEOUT_S)
                _write_status(run, state="done", finished=time.time(), cpu_s=out["summary"]["cpu_s"])
        except Exception as e:                  # noqa: BLE001
            _write_status(run, state="error", finished=time.time(), error=str(e)[-1500:])
    POOL.submit(job)
    return rid


def list_runs(ws: Workspace) -> list[dict]:
    out = []
    for d in sorted(ws.runs.iterdir(), reverse=True):
        if d.is_dir():
            st = _status(d)
            out.append({"id": d.name, **{k: st.get(k) for k in ("state", "name", "created", "cpu_s", "progress",
                                                                 "error", "model")}})
    return out


def run_results(ws: Workspace, rid: str) -> dict:
    from .run import load_run
    run = ws.runs / _safe_id(rid)
    if run.is_symlink() or not run.is_dir():
        raise KeyError(rid)
    r = load_run(run)
    res, g = r["results"], r["grid"]
    return {"id": rid, "name": r["case"]["case"]["name"], "summary": r["summary"],
            "time_s": res["time_s"].tolist(), "met": _round(res["met"]),
            "bounds_um": g["bounds_um"], "class_labels": g["class_labels"],
            "species": g["species"]["aerosol"], "gases": g["species"]["gases"],
            "number_m3": _round(res["number_m3"]), "mass_ugm3": _round(res["mass_ugm3"]),
            "gas_ugm3": _round(res["gas_ugm3"]), "post": r.get("post"),
            "het": _round(res["het"]) if "het" in res else None,
            "kernel": r["case"]["kernel"]["isorropia"], "notes": g.get("notes", []),
            "model": "partmc" if (run / "pmc.spec").exists() else "scram",
            "truth": json.loads((run / "truth.json").read_text()) if (run / "truth.json").exists() else None,
            "partmc_options": (json.loads((run / "partmc_options.json").read_text())
                               if (run / "partmc_options.json").exists() else None)}


class RunNotReady(RuntimeError):
    pass


class ExportTooLarge(ValueError):
    pass


def _export_regular_file(name: str, *, directory_fd: int) -> bytes:
    fd = os.open(name, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK, dir_fd=directory_fd)
    with os.fdopen(fd, "rb") as source:
        info = os.fstat(source.fileno())
        if not stat.S_ISREG(info.st_mode):
            raise ValueError("export source is not a regular file")
        if info.st_size > MAX_EXPORT_BODY:
            raise ExportTooLarge("raw result file exceeds 128 MiB")
        data = source.read(MAX_EXPORT_BODY + 1)
    if len(data) > MAX_EXPORT_BODY:
        raise ExportTooLarge("raw result file exceeds 128 MiB")
    return data


def run_export(ws: Workspace, rid: str) -> bytes:
    """Export a completed owned SCRAM or PartMC run with bounded raw truth."""
    from .export_bundle import result_archive
    run = ws.runs / _safe_id(rid)
    if run.is_symlink() or not run.is_dir():
        raise KeyError(rid)
    if _status(run).get("state") != "done":
        raise RunNotReady("run must finish before export")
    directory_fd = os.open(run, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW)
    raw_files: dict[str, bytes] = {}
    try:
        case = json.loads(_export_regular_file("case.json", directory_fd=directory_fd))
        for name in ("results.npz", "grid.json", "case.resolved.json", "post.json",
                     "truth.json", "partmc_options.json"):
            if name in {"post.json", "truth.json", "partmc_options.json"} and not (run / name).exists():
                continue
            raw_files[name] = _export_regular_file(name, directory_fd=directory_fd)
        out_path = run / "out"
        if out_path.exists():
            out_fd = os.open("out", os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW, dir_fd=directory_fd)
            try:
                names = sorted(name for name in os.listdir(out_fd)
                               if re.fullmatch(r"pmc_0001_[0-9]+\.nc", name))
                if len(names) > MAX_EXPORT_FILES:
                    raise ExportTooLarge("too many PartMC NetCDF files for one export")
                for name in names:
                    raw_files["out/" + name] = _export_regular_file(name, directory_fd=out_fd)
            finally:
                os.close(out_fd)
        if sum(len(data) for data in raw_files.values()) > MAX_EXPORT_BODY:
            raise ExportTooLarge("raw result files exceed 128 MiB")
    finally:
        os.close(directory_fd)
    try:
        archive = result_archive(rid, case, run_results(ws, rid), raw_files)
    except ValueError as exc:
        if "too large" in str(exc):
            raise ExportTooLarge(str(exc)) from exc
        raise
    if len(archive) > MAX_EXPORT_BODY:
        raise ExportTooLarge("export archive exceeds 128 MiB")
    return archive


def run_post(ws: Workspace, rid: str, cfg: dict) -> dict:
    from .post import run_all
    from .run import load_run
    run = ws.runs / _safe_id(rid)
    if run.is_symlink() or not run.is_dir():
        raise KeyError(rid)
    r = load_run(run)
    return run_all(r, cfg, run)


# ------------------------------------------------------------------ experiments
def start_experiment(ws: Workspace, spec: dict) -> str:
    from .experiments import expand, run_set
    if shared_mode() and len(spec.get("variants", [])) > MAX_VARIANTS:
        raise ValueError(f"at most {MAX_VARIANTS} variants per experiment set")
    from .experiments import variant_model
    for v, c in zip(spec.get("variants", []), expand(spec)):
        model, opts = variant_model(v)
        if model == "partmc":
            partmc_options(c, opts)
        else:
            check_cost(c)
    _check_quota(ws)
    eid = _new_id(spec.get("name", "experiment"))
    d = ws.exps / eid
    d.mkdir()
    (d / "status.json").write_text(json.dumps({"state": "queued", "variants": []}))
    POOL.submit(lambda: run_set(spec, d, workers=int(os.environ.get("SCRAMBOX_EXP_WORKERS", "3"))))
    return eid


def list_experiments(ws: Workspace) -> list[dict]:
    out = []
    for d in sorted(ws.exps.iterdir(), reverse=True):
        if not d.is_dir():
            continue
        st = json.loads((d / "status.json").read_text()) if (d / "status.json").exists() else {}
        spec = json.loads((d / "experiment.json").read_text()) if (d / "experiment.json").exists() else {}
        out.append({"id": d.name, "name": spec.get("name"), "state": st.get("state"),
                    "variants": [v.get("label") for v in spec.get("variants", [])]})
    return out


def experiment(ws: Workspace, eid: str) -> dict:
    from .experiments import compare
    d = ws.exps / _safe_id(eid)
    if not d.is_dir():
        raise KeyError(eid)
    st = json.loads((d / "status.json").read_text())
    if st.get("state") in ("queued", "running") or not (d / "experiment.json").exists():
        return {"id": eid, "status": st.get("state"), "progress": st.get("variants", [])}
    out = compare(d)
    out["id"] = eid
    return out


# ------------------------------------------------------------------ trusted website proxy
def trusted_proxy_secret() -> bytes | None:
    path = Path(os.environ.get("SCRAMBOX_TRUSTED_PROXY_SECRET_FILE",
                               data_dir() / "trusted_proxy_secret"))
    if not path.is_file():
        return None
    secret = path.read_bytes().strip()
    if len(secret) < 32:
        raise RuntimeError("trusted proxy secret must contain at least 32 bytes")
    return secret


# ------------------------------------------------------------------ HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = "scrambox/0.3"
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):         # quiet
        pass

    def _https(self) -> bool:
        return self.headers.get("X-Forwarded-Proto", "").lower() == "https" or \
            "Tailscale-Funnel-Request" in self.headers

    def _send(self, code, body, ctype="application/json; charset=utf-8", cookie=None, download_name=None):
        data = body if isinstance(body, bytes) else json.dumps(body, ensure_ascii=False, allow_nan=False,
                                                             default=_json_default).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("X-Frame-Options", "SAMEORIGIN")
        if cookie is not None:
            self.send_header("Set-Cookie", cookie)
        if download_name is not None:
            self.send_header("Content-Disposition", f'attachment; filename="{download_name}"')
        self.end_headers()
        self.wfile.write(data)

    def _cookie(self, token: str | None) -> str:
        attrs = "Path=/; HttpOnly; SameSite=Lax" + ("; Secure" if self._https() else "")
        if token is None:
            return f"sb_session=; Max-Age=0; {attrs}"
        return f"sb_session={token}; Max-Age={SESSION_DAYS * 86400}; {attrs}"

    def _token(self) -> str | None:
        c = SimpleCookie(self.headers.get("Cookie", ""))
        return c["sb_session"].value if "sb_session" in c else None

    def _ip(self) -> str:
        return (self.headers.get("X-Forwarded-For", "") or self.client_address[0]).split(",")[0].strip()

    def _body_bytes(self) -> bytes:
        try:
            n = int(self.headers.get("Content-Length") or 0)
        except ValueError:
            raise ValueError("invalid content length")
        if n < 0 or n > MAX_BODY:
            raise ValueError("request too large")
        return self.rfile.read(n)

    def _trusted_proxy_user(self, method: str, raw_body: bytes) -> str | None:
        user = self.headers.get("X-Scram-User")
        timestamp = self.headers.get("X-Scram-Timestamp")
        signature = self.headers.get("X-Scram-Signature")
        if not any((user, timestamp, signature)):
            return None
        if not all((user, timestamp, signature)):
            raise PermissionError("incomplete trusted proxy signature")
        if self._token():
            raise PermissionError("trusted proxy request cannot use session cookie")
        if not re.fullmatch(r"5002:[1-9][0-9]*", user):
            raise PermissionError("invalid trusted proxy user")
        if not re.fullmatch(r"[0-9]{10,12}", timestamp):
            raise PermissionError("invalid trusted proxy timestamp")
        if abs(time.time() - int(timestamp)) > 60:
            raise PermissionError("expired trusted proxy signature")
        if not re.fullmatch(r"[0-9a-f]{64}", signature):
            raise PermissionError("invalid trusted proxy signature")
        secret = trusted_proxy_secret()
        if secret is None:
            raise PermissionError("trusted proxy is not configured")
        message = "\n".join((user, timestamp, method.upper(), self.path,
                             hashlib.sha256(raw_body).hexdigest())).encode()
        if not hmac.compare_digest(hmac.new(secret, message, hashlib.sha256).hexdigest(), signature):
            raise PermissionError("invalid trusted proxy signature")
        return user

    def _route(self, method):
        path = unquote(urlparse(self.path).path)
        parts = [p for p in path.split("/") if p]
        try:
            if not parts or parts[0] != "api":
                return self._static(path)
            p = parts[1:]
            raw_body = self._body_bytes()
            trusted_user = self._trusted_proxy_user(method, raw_body)
            body = json.loads(raw_body or b"{}") if method == "POST" else None
            shared = shared_mode()
            user = trusted_user if trusted_user is not None else (session_user(self._token()) if shared else None)
            if p == ["session"]:
                mode = "trusted" if trusted_user is not None else ("shared" if shared else "local")
                return self._send(200, {"mode": mode, "user": user})
            if p == ["login"] and method == "POST":
                if trusted_user is not None:
                    return self._send(403, {"error": "website identity is managed by the website"})
                if not shared:
                    return self._send(200, {"mode": "local", "user": None})
                token = login(body.get("code", ""), body.get("name", ""), self._ip())
                return self._send(200, {"mode": "shared", "user": body.get("name", "").strip()}, cookie=self._cookie(token))
            if p == ["logout"] and method == "POST":
                if trusted_user is not None:
                    return self._send(403, {"error": "website identity is managed by the website"})
                logout(self._token())
                return self._send(200, {"ok": True}, cookie=self._cookie(None))
            if shared and user is None:
                return self._send(401, {"error": "sign in required", "login": True})
            ws = workspace_for(user)
            base = ws.root
            if method == "GET":
                if p == ["meta"]:
                    return self._send(200, meta())
                if p == ["scenarios"]:
                    return self._send(200, scn.index(base))
                if len(p) == 2 and p[0] == "scenarios":
                    return self._send(200, scn.get(p[1], base))
                if p == ["runs"]:
                    return self._send(200, list_runs(ws))
                if len(p) == 2 and p[0] == "runs":
                    return self._send(200, _status(ws.runs / _safe_id(p[1])))
                if len(p) == 3 and p[0] == "runs" and p[2] == "results":
                    return self._send(200, run_results(ws, p[1]))
                if len(p) == 3 and p[0] == "runs" and p[2] == "export":
                    rid = _safe_id(p[1])
                    model = _status(ws.runs / rid).get("model")
                    prefix = "partmcbox" if model == "partmc" else "scrambox"
                    return self._send(200, run_export(ws, rid), "application/zip",
                                      download_name=f"{prefix}-{rid}.zip")
                if len(p) == 3 and p[0] == "runs" and p[2] == "particles":
                    from urllib.parse import parse_qs
                    from .partmc import particles_at
                    q = parse_qs(urlparse(self.path).query)
                    particle_run = ws.runs / _safe_id(p[1])
                    if particle_run.is_symlink() or not particle_run.is_dir():
                        raise KeyError(p[1])
                    return self._send(200, particles_at(particle_run, int(q.get("k", ["-1"])[0])
                                                        if q.get("k", ["-1"])[0] != "-1" else 10 ** 6))
                if len(p) == 3 and p[0] == "runs" and p[2] == "case":
                    f = ws.runs / _safe_id(p[1]) / "case.json"
                    if not f.exists():
                        raise KeyError(p[1])
                    return self._send(200, json.loads(f.read_text()))
                if p == ["experiments"]:
                    return self._send(200, list_experiments(ws))
                if len(p) == 2 and p[0] == "experiments":
                    return self._send(200, experiment(ws, p[1]))
            elif method == "POST":
                if p == ["describe"]:
                    return self._send(200, scn.describe(body))
                if p == ["grid"]:
                    return self._send(200, grid_info(body))
                if p == ["resolve"]:
                    return self._send(200, casemod.resolve(body))
                if p == ["scenarios"]:
                    return self._send(200, scn.save_user(body, base))
                if p == ["runs"]:
                    return self._send(200, {"id": start_run(ws, body["case"], body.get("name"),
                                                            body.get("model", "scram"), body.get("partmc"))})
                if len(p) == 3 and p[0] == "runs" and p[2] == "post":
                    return self._send(200, run_post(ws, p[1], body))
                if p == ["experiments"]:
                    return self._send(200, {"id": start_experiment(ws, body)})
            elif method == "DELETE":
                if len(p) == 2 and p[0] == "scenarios":
                    scn.delete_user(p[1], base)
                    return self._send(200, {"deleted": p[1]})
                if len(p) == 2 and p[0] in ("runs", "experiments"):
                    d = (ws.runs if p[0] == "runs" else ws.exps) / _safe_id(p[1])
                    if d.is_symlink() or not d.is_dir():
                        raise KeyError(p[1])
                    shutil.rmtree(d)
                    return self._send(200, {"deleted": p[1]})
            return self._send(404, {"error": "not found"})
        except PermissionError as e:
            return self._send(403, {"error": str(e)})
        except RunNotReady as e:
            return self._send(409, {"error": str(e)})
        except ExportTooLarge as e:
            return self._send(413, {"error": str(e)})
        except KeyError as e:
            return self._send(404, {"error": f"not found: {e}"})
        except (ValueError, FileNotFoundError, json.JSONDecodeError) as e:
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
    httpd.daemon_threads = True
    print(f"SCRAM box app on http://{host}:{port}  mode={'shared' if shared_mode() else 'local'}  "
          f"data={data_dir()}", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        pass
