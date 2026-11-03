"""Dust follow-up: does the outer coupling step (dt_s, RH/T update interval) or the wet-diameter rule explain the lagging growth?"""
import copy, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.case import set_path
from scrambox.post import run_all
from scrambox.run import load_run, run_case
OUT = ROOT / "runs" / "dust_diag2"; OUT.mkdir(parents=True, exist_ok=True)
POST = {"ccn": {"enabled": True, "every": 2}}
V = {"x8_dt600": {}, "x8_dt60": {"time.dt_s": 60.0}, "x8_dt60_sub1": {"time.dt_s": 60.0, "time.substep_s": 1.0},
     "x8_dt600_gerber": {"processes.wet_diameter": "gerber"}, "x8_dt60_nocoag": {"time.dt_s": 60.0, "processes.coagulation": False},
     "x1_dt60": {"time.dt_s": 60.0, "_k": 1}, "x8_dt60_fixeddens": {"time.dt_s": 60.0, "processes.density.mode": "fixed"},
     "x8_dt600_sub1": {"time.substep_s": 1.0}, "x8_dt600_sub60": {"time.substep_s": 60.0}}
def go(name):
    ov = dict(V[name]); k = ov.pop("_k", 8)
    c = copy.deepcopy(scn.get("dust_crustal")["case"]); c["mixing"] = {"preset": "NC1"}; c["postprocess"] = POST
    c["size_bins"] = {"mode": "log", "n": 14 * k, "dmin_um": 0.01, "dmax_um": 10.0}
    for p, v in ov.items(): c = set_path(c, p, v)
    d = OUT / name
    try:
        run_case(c, d, timeout_s=40000); run_all(load_run(d), POST, d); r = load_run(d)
        g = {str(x["S_pct"]): round(x["ccn_cm3"], 1) for x in r["post"]["ccn"]["resolved"][-1]["spectrum"]}
        sp = r["grid"]["species"]["aerosol"]; m = r["results"]["mass_ugm3"][-1].sum(axis=(0, 1))
        return name, {"ccn0.1": g.get("0.1"), "ccn0.2": g.get("0.2"), "NO3": round(float(m[sp.index("PNO3")]), 3), "cpu_s": round(r["summary"].get("cpu_s") or 0, 1)}
    except Exception as e:
        return name, {"error": str(e)[-300:]}
res = {}
with ThreadPoolExecutor(9) as ex:
    for name, v in ex.map(go, V):
        res[name] = v; print(name, v, flush=True); (OUT / "diag2.json").write_text(json.dumps(res, indent=1))
print("=== dust diag2 done (truth seed1 CCN0.1 1487.1; diag truth 1494.4)")
