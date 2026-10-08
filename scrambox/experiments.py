"""Experiment sets: one base case and several variants (path overrides), run together
and compared on the same metrics, e.g. NC1 vs NC6 vs NC11, or ISORROPIA v1 vs v2.3.

spec = {"name", "base": <case>, "variants": [{"label", "overrides": {"mixing.preset": "NC1"}}],
        "postprocess": <case.postprocess or None>}
"""
from __future__ import annotations

import copy
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

from . import case as casemod

PRESET_DIMENSIONS = {
    "mixing.preset": ["NC1", "NC6", "NC11"],
    "kernel.isorropia": ["none", "v1", "v23"],
    "processes.thermodynamics": ["dynamic", "equilibrium"],
    "time.substep_s": [600.0, 60.0, 10.0],
    "processes.coagulation": [True, False],
}


def expand(spec: dict) -> list[dict]:
    cases = []
    for v in spec["variants"]:
        c = copy.deepcopy(spec["base"])
        for path, val in v.get("overrides", {}).items():
            c = casemod.set_path(c, path, val)
        if spec.get("postprocess") is not None:
            c["postprocess"] = copy.deepcopy(spec["postprocess"])
        c.setdefault("case", {})["name"] = f"{spec.get('name', 'experiment')} - {v['label']}"
        cases.append(c)
    return cases


def _slug(s):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s).strip("_")[:60] or "variant"


def run_set(spec: dict, exp_dir: str | Path, workers: int = 4, on_update=None) -> dict:
    from .run import run_case
    exp = Path(exp_dir)
    exp.mkdir(parents=True, exist_ok=True)
    (exp / "experiment.json").write_text(json.dumps(spec, indent=1, ensure_ascii=False))
    cases = expand(spec)
    status = {"state": "running", "variants": [{"label": v["label"], "state": "queued"} for v in spec["variants"]],
              "started": time.time()}

    def save():
        (exp / "status.json").write_text(json.dumps(status, indent=1))
        if on_update:
            on_update(status)

    save()

    def one(i):
        status["variants"][i]["state"] = "running"
        save()
        d = exp / f"{i:02d}_{_slug(spec['variants'][i]['label'])}"
        try:
            out = run_case(cases[i], d)
            status["variants"][i].update(state="done", run_dir=str(d), cpu_s=out["summary"]["cpu_s"])
        except Exception as e:                      # noqa: BLE001 - reported to the user
            status["variants"][i].update(state="error", error=str(e)[-800:], run_dir=str(d))
        save()

    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        list(ex.map(one, range(len(cases))))
    status["state"] = "done" if all(v["state"] == "done" for v in status["variants"]) else "error"
    status["finished"] = time.time()
    save()
    return status


def metrics(run: dict) -> dict:
    """Time series used to compare variants."""
    r, g = run["results"], run["grid"]
    sp = g["species"]["aerosol"]
    iw = sp.index("PH2O")
    m = r["mass_ugm3"]                                   # (t, b, c, s)
    upper = np.array(g["bounds_um"][1:])
    dry = m.sum(axis=(1, 2))                             # (t, s)
    pm1 = m[:, upper <= 1.0 + 1e-9].sum(axis=(1, 2))     # (t, s)
    out = {"time_h": (r["time_s"] / 3600.0).tolist(),
           "N_cm3": (r["number_m3"].sum(axis=(1, 2)) / 1e6).tolist(),
           "PM_dry": (dry.sum(axis=1) - dry[:, iw]).tolist(),
           "PM1_dry": (pm1.sum(axis=1) - pm1[:, iw]).tolist(),
           "water": dry[:, iw].tolist(),
           "species": {s: dry[:, i].tolist() for i, s in enumerate(sp) if dry[:, i].max() > 0}}
    post = run.get("post") or {}
    if "optics" in post:
        o = post["optics"]
        out["post_time_h"] = [t / 3600.0 for t in o["time_s"]]
        out["ext550"] = [x["ext550_Mm"] for x in o["resolved"]]
        out["ext550_internal"] = [x["ext550_Mm"] for x in o["internal_mean"]]
        out["ssa550"] = [x["ssa550"] for x in o["resolved"]]
        out["bc_abs_enh"] = [x["bc_abs_enhancement"] for x in o["resolved"]]
    if "ccn" in post:
        c = post["ccn"]
        out["post_time_h"] = [t / 3600.0 for t in c["time_s"]]
        out["ccn"] = {str(s["S_pct"]): [st["spectrum"][k]["ccn_cm3"] for st in c["resolved"]]
                      for k, s in enumerate(c["resolved"][0]["spectrum"])}
        out["nd"] = {str(a["w_m_s"]): [st["activation"][k]["nd_cm3"] for st in c["resolved"]]
                     for k, a in enumerate(c["resolved"][0]["activation"])}
    if "hetchem" in post:
        h = post["hetchem"]
        out["post_time_h"] = [t / 3600.0 for t in h["time_s"]]
        out["k_het"] = {}
        for gas, by in h["states"][0]["rules"].items():
            for rule in by["resolved"]:
                out["k_het"][f"{gas}:{rule}"] = [st["rules"][gas]["resolved"][rule]["k_s1"] for st in h["states"]]
                out["k_het"][f"{gas}:{rule}:internal"] = [st["rules"][gas]["internal_mean"][rule]["k_s1"]
                                                         for st in h["states"]]
    return out


def compare(exp_dir: str | Path) -> dict:
    from .run import load_run
    exp = Path(exp_dir)
    spec = json.loads((exp / "experiment.json").read_text())
    status = json.loads((exp / "status.json").read_text())
    rows = []
    for v in status["variants"]:
        row = {"label": v["label"], "state": v["state"]}
        if v["state"] == "done":
            row["metrics"] = metrics(load_run(v["run_dir"]))
        elif "error" in v:
            row["error"] = v["error"]
        rows.append(row)
    return {"name": spec.get("name"), "status": status["state"], "variants": rows}
