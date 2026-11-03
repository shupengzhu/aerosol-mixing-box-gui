"""Dust follow-up 3: redistribution scheme x resolution, and convergence to x16/x32 (one class)."""
import copy, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.case import set_path
from scrambox.post import run_all
from scrambox.run import load_run, run_case
OUT = ROOT / "runs" / "dust_diag3"; OUT.mkdir(parents=True, exist_ok=True)
POST = {"ccn": {"enabled": True, "every": 2}}
V = {}
for red in ("moving_center_dualpivot", "core_fluxform", "core_ppm_cutnum", "core_conserv", "legacy", "core_smallgrow"):
    for k in (1, 8):
        V[f"{red}_x{k}"] = {"processes.redistribution": red, "_k": k}
V["moving_center_dualpivot_x16"] = {"_k": 16}; V["moving_center_dualpivot_x32"] = {"_k": 32}
V["core_fluxform_x16"] = {"processes.redistribution": "core_fluxform", "_k": 16}
V["core_fluxform_x8_nocoag"] = {"processes.redistribution": "core_fluxform", "_k": 8, "processes.coagulation": False}
def go(name):
    ov = dict(V[name]); k = ov.pop("_k", 8)
    c = copy.deepcopy(scn.get("dust_crustal")["case"]); c["mixing"] = {"preset": "NC1"}; c["postprocess"] = POST
    c["size_bins"] = {"mode": "log", "n": 14 * k, "dmin_um": 0.01, "dmax_um": 10.0}
    for p, v in ov.items(): c = set_path(c, p, v)
    d = OUT / name
    try:
        run_case(c, d, timeout_s=80000); run_all(load_run(d), POST, d); r = load_run(d)
        g = {str(x["S_pct"]): round(x["ccn_cm3"], 1) for x in r["post"]["ccn"]["resolved"][-1]["spectrum"]}
        sp = r["grid"]["species"]["aerosol"]; m = r["results"]["mass_ugm3"][-1].sum(axis=(0, 1))
        n = r["results"]["number_m3"][-1].sum(1) / 1e6; b = np.array(r["grid"]["bounds_um"]); mid = np.sqrt(b[1:] * b[:-1])
        return name, {"N": round(float(n.sum()), 1), "ccn0.1": g.get("0.1"), "ccn0.2": g.get("0.2"), "NO3": round(float(m[sp.index("PNO3")]), 3),
                      "N_0.04_0.107": round(float(n[(mid > 0.04) & (mid < 0.107)].sum()), 1), "N_0.107_0.3": round(float(n[(mid >= 0.107) & (mid < 0.3)].sum()), 1),
                      "cpu_s": round(r["summary"].get("cpu_s") or 0, 1)}
    except Exception as e:
        return name, {"error": str(e)[-300:]}
res = {}
with ThreadPoolExecutor(10) as ex:
    for name, v in ex.map(go, V):
        res[name] = v; print(name, v, flush=True); (OUT / "diag3.json").write_text(json.dumps(res, indent=1))
print("=== dust diag3 done (PartMC x8 projection: N_0.04_0.107 = 0.0, N_0.107_0.3 = 1881.6, CCN0.1 1494)")
