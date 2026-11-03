"""Dust follow-up 4: semi-volatile nitrate split between two narrow sulfate modes (0.06 and 0.15 um), rising RH,
no dust, no coagulation. SCRAM one class at x8 vs PartMC 5000 particles. Isolates size-dependent HNO3/NH3 uptake."""
import copy, json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.case import set_path
from scrambox.partmc import run_partmc, project
from scrambox.post import run_all
from scrambox.run import load_run, run_case
OUT = ROOT / "runs" / "dust_diag4"; OUT.mkdir(parents=True, exist_ok=True)
c = copy.deepcopy(scn.get("dust_crustal")["case"]); c["mixing"] = {"preset": "NC1"}; c["postprocess"] = {"ccn": {"enabled": True, "every": 2}}
c["initial"]["particles"]["modes"] = [{"name": "small sulfate", "number_cm3": 1000, "dg_um": 0.06, "sigma_g": 1.08, "mass_fraction": {"PSO4": 0.727, "PNH4": 0.273}},
                                      {"name": "large sulfate", "number_cm3": 1000, "dg_um": 0.15, "sigma_g": 1.08, "mass_fraction": {"PSO4": 0.727, "PNH4": 0.273}}]
c = set_path(c, "processes.coagulation", False)
c["size_bins"] = {"mode": "log", "n": 14 * 8, "dmin_um": 0.01, "dmax_um": 10.0}
def binned(d):
    r = load_run(d); g = r["grid"]; b = np.array(g["bounds_um"]); mid = np.sqrt(b[1:] * b[:-1]); sp = g["species"]["aerosol"]
    out = {"t_h": (r["results"]["time_s"] / 3600).tolist(), "small": [], "large": []}
    for k in range(len(r["results"]["time_s"])):
        n = r["results"]["number_m3"][k].sum(1) / 1e6; m = r["results"]["mass_ugm3"][k].sum(1)
        for lab, sel in (("small", mid < 0.1), ("large", mid >= 0.1)):
            dry = m[sel].sum() - m[sel][:, sp.index("PH2O")].sum()
            dmean = float((n[sel] * mid[sel]).sum() / max(n[sel].sum(), 1e-9))
            out[lab].append({"N": float(n[sel].sum()), "NO3": float(m[sel][:, sp.index("PNO3")].sum()), "NH4": float(m[sel][:, sp.index("PNH4")].sum()),
                             "H2O": float(m[sel][:, sp.index("PH2O")].sum()), "dry": float(dry), "D_mean": dmean})
    cc = r["post"]["ccn"]["resolved"][-1]["spectrum"]; out["ccn_end"] = {str(x["S_pct"]): x["ccn_cm3"] for x in cc}
    return out
res = {}
d = OUT / "scram_x8"; run_case(c, d, timeout_s=40000); run_all(load_run(d), c["postprocess"], d); res["scram_x8"] = binned(d)
c1 = copy.deepcopy(c); c1["size_bins"]["n"] = 14 * 32; d = OUT / "scram_x32"; run_case(c1, d, timeout_s=80000); run_all(load_run(d), c["postprocess"], d); res["scram_x32"] = binned(d)
pd_ = OUT / "pmc"; cp = copy.deepcopy(c); cp["size_bins"]["n"] = 14; run_partmc(cp, pd_, n_part=5000, seed=1)
pp = OUT / "pmc_x8"; project(pd_, pp, {"preset": "NC1"}, {"size_bins": c["size_bins"]}); run_all(load_run(pp), c["postprocess"], pp); res["pmc_x8"] = binned(pp)
(OUT / "diag4.json").write_text(json.dumps(res, default=float))
for k, v in res.items():
    print(k, "end small N/NO3/D:", round(v["small"][-1]["N"], 1), round(v["small"][-1]["NO3"], 3), round(v["small"][-1]["D_mean"], 4),
          "large N/NO3/D:", round(v["large"][-1]["N"], 1), round(v["large"][-1]["NO3"], 3), round(v["large"][-1]["D_mean"], 4), "CCN0.1", round(v["ccn_end"]["0.1"], 1))
print("=== diag4 done")
