"""E1: urban plume, NC6 vs NC9BC (explicit BC classes) at x1/x2/x4/x8 size resolution:
SCRAM simulation and PartMC (seed 1) projected on the same grid; CCN, E_abs, chi vs particle-resolved truth."""
import copy, json, sys, time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.partmc import chi_sections, project
from scrambox.post import run_all
from scrambox.run import load_run, run_case
SC = "urban_plume_aging"; PM = ROOT / "runs" / "wp2" / SC / "partmc_s1"
OUT = ROOT / "runs" / "e1_nc9bc_binres"; OUT.mkdir(parents=True, exist_ok=True)
POST = {"optics": {"enabled": True, "every": 8}, "ccn": {"enabled": True, "every": 8}}
def bins(n): return {"mode": "log", "n": n, "dmin_um": 0.005, "dmax_um": 2.5}
def summ(d):
    r = load_run(d); m, n = r["results"]["mass_ugm3"][-1], r["results"]["number_m3"][-1]; sp = r["grid"]["species"]["aerosol"]
    o = r["post"]["optics"]["resolved"][-1]; cc = {str(x["S_pct"]): round(x["ccn_cm3"], 1) for x in r["post"]["ccn"]["resolved"][-1]["spectrum"]}
    return {"n_class": len(r["grid"]["class_labels"]), "nsize": len(r["grid"]["bounds_um"]) - 1, "chi": round(chi_sections(m, n, sp), 4),
            "eabs": round(o["bc_abs_enhancement"], 3), "ccn0.1": cc.get("0.1"), "ccn0.2": cc.get("0.2"), "ccn0.5": cc.get("0.5"),
            "cpu_s": round((r.get("summary") or {}).get("cpu_s") or 0, 1), "kernel": (r.get("summary") or {}).get("kernel")}
def scram(mix, k):
    c = copy.deepcopy(scn.get(SC)["case"]); c["mixing"] = {"preset": mix}; c["size_bins"] = bins(16 * k); c["postprocess"] = POST
    d = OUT / f"scram_{mix}_x{k}"; run_case(c, d, timeout_s=40000); run_all(load_run(d), POST, d); return f"scram_{mix}_x{k}", summ(d)
def pmc(mix, k):
    d = OUT / f"pmc_{mix}_x{k}"; project(PM, d, {"preset": mix}, {"size_bins": bins(16 * k)}); run_all(load_run(d), POST, d); return f"pmc_{mix}_x{k}", summ(d)
res = {}
t = json.loads((PM / "truth.json").read_text())
res["truth"] = {"chi": round(t["chi"][-1], 4), **{f"ccn{s}": round(t["ccn_cm3"][s][-1], 1) for s in ("0.1", "0.2", "0.5")}}
jobs = [(f, mix, k) for mix in ("NC6", "NC9BC") for k in (1, 2, 4, 8) for f in (scram, pmc)]
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as ex:
    futs = [ex.submit(f, mix, k) for f, mix, k in jobs]
    for fu in futs:
        try:
            k, v = fu.result(); res[k] = v; print(k, v, flush=True)
        except Exception as e:
            print("ERROR", str(e)[-300:], flush=True)
        (OUT / "e1.json").write_text(json.dumps(res, indent=1))
print("=== e1 done")
