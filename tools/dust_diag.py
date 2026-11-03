"""Dust CCN(0.1%) deficit: process / composition / size / time decomposition (SCRAM one class vs PartMC seed 1).
Each variant runs SCRAM NC1 at x1 and x8 size resolution and a matching PartMC run (5000 particles) with the same switches.
Usage: python tools/dust_diag.py <workers>   -> runs/dust_diag/dust_diag.json"""
import copy, json, sys, os
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.case import set_path
from scrambox.partmc import run_partmc, truth, project, read_particles
from scrambox.post import run_all
from scrambox.run import load_run, run_case
OUT = ROOT / "runs" / "dust_diag"; OUT.mkdir(parents=True, exist_ok=True)
base = copy.deepcopy(scn.get("dust_crustal")["case"])
POST = {"ccn": {"enabled": True, "every": 1}}
def variant(name):
    c = copy.deepcopy(base); c["mixing"] = {"preset": "NC1"}; c["postprocess"] = POST
    if name == "full": pass
    elif name == "nocoag": c = set_path(c, "processes.coagulation", False)
    elif name == "nocond": c = set_path(c, "processes.condensation", False)
    elif name == "no_hno3": c["initial"]["gas_ugm3"].pop("HNO3", None)
    elif name == "no_nh3": c["initial"]["gas_ugm3"].pop("NH3", None)
    elif name == "no_hcl": c["initial"]["gas_ugm3"].pop("HCL", None)
    elif name == "no_coarse_dust": c["initial"]["particles"]["modes"] = [m for m in c["initial"]["particles"]["modes"] if m["name"] != "dust coarse"]
    elif name == "no_dust": c["initial"]["particles"]["modes"] = [m for m in c["initial"]["particles"]["modes"] if not m["name"].startswith("dust")]
    elif name == "no_sulf_src": c["emissions"]["gas_ugm3_s"] = {}
    elif name == "v1_inert_dust": c["kernel"] = {"isorropia": "v1"}
    elif name == "const_rh": c["environment"]["time_series"]["relative_humidity"] = [0.85] * len(c["environment"]["time_series"]["time_s"]); c["environment"]["relative_humidity"] = 0.85
    elif name == "fixed_density": c = set_path(c, "processes.density.mode", "fixed")
    elif name == "wet_kappa": c = set_path(c, "processes.wet_diameter", "kappa")
    else: raise ValueError(name)
    return c
NAMES = ["full", "nocoag", "nocond", "no_hno3", "no_nh3", "no_hcl", "no_coarse_dust", "no_dust", "no_sulf_src", "v1_inert_dust", "const_rh"]
def ccn_series(d):
    r = load_run(d); c = r["post"]["ccn"]
    return {"t_h": [t / 3600 for t in c["time_s"]], **{str(s["S_pct"]): [st["spectrum"][k]["ccn_cm3"] for st in c["resolved"]] for k, s in enumerate(c["resolved"][0]["spectrum"]) if s["S_pct"] in (0.1, 0.2, 0.5)}}
def size_kappa(d):
    """end-of-run per-bin number (cm-3) and volume-weighted kappa (one class)."""
    from scrambox.post.common import populations, KAPPA
    r = load_run(d); pop = populations(r, len(r["results"]["time_s"]) - 1)
    kap = np.array([KAPPA.get(s, 0.1) for s in pop.species]); m = np.maximum(pop.mass, 0); m[..., pop.iw] = 0
    vol = m / pop.density_gcm3; vt = vol.sum(-1)
    kk = (vol * kap).sum(-1) / np.where(vt > 0, vt, 1)
    return {"bounds_um": list(map(float, pop.bounds_um)), "n_cm3": (pop.number.sum(1) / 1e6).tolist(), "kappa": kk[:, 0].tolist() if pop.nc == 1 else (kk * pop.number).sum(1).tolist()}
def go(name):
    out = {}
    try:
        for k in (1, 8):
            c = variant(name); c["size_bins"] = {"mode": "log", "n": 14 * k, "dmin_um": 0.01, "dmax_um": 10.0}
            d = OUT / f"{name}_x{k}"; run_case(c, d, timeout_s=40000); run_all(load_run(d), POST, d)
            out[f"scram_x{k}"] = ccn_series(d); out[f"scram_x{k}_size"] = size_kappa(d)
        c = variant(name); c["size_bins"] = {"mode": "log", "n": 14, "dmin_um": 0.01, "dmax_um": 10.0}
        if c["kernel"]["isorropia"] == "v1":
            c["kernel"]["isorropia"] = "v1"
        pd_ = OUT / f"{name}_pmc"; run_partmc(c, pd_, n_part=5000, seed=1)
        t = truth(pd_); out["truth"] = {"t_h": [x / 3600 for x in t["time_s"]], **{s: t["ccn_cm3"][s] for s in ("0.1", "0.2", "0.5")}}
        for k in (1, 8):
            pp = OUT / f"{name}_pmc_x{k}"; project(pd_, pp, {"preset": "NC1"}, {"size_bins": {"mode": "log", "n": 14 * k, "dmin_um": 0.01, "dmax_um": 10.0}}); run_all(load_run(pp), POST, pp)
            out[f"pmc_x{k}"] = ccn_series(pp); out[f"pmc_x{k}_size"] = size_kappa(pp)
        print(name, "end ratios x1/x8/truth:", round(out["scram_x1"]["0.1"][-1], 1), round(out["scram_x8"]["0.1"][-1], 1), round(out["truth"]["0.1"][-1], 1), flush=True)
    except Exception as e:
        out["error"] = str(e)[-400:]; print(name, "ERROR", out["error"], flush=True)
    return name, out
res = {}
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as ex:
    for name, v in ex.map(go, NAMES):
        res[name] = v; (OUT / "dust_diag.json").write_text(json.dumps(res, default=float))
print("=== dust diag done")
