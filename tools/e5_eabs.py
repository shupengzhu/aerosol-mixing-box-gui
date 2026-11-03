"""E_abs consistency (urban, 24 h): every class configuration evaluated with BOTH optics paths
(default post-processing and prognostic-moment core-shell LUT used for the particle-resolved truth)."""
import copy, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn
from scrambox.partmc import chi_sections, project, truth_optics
from scrambox.post import optics, run_all
from scrambox.post.common import populations
from scrambox.run import load_run, run_case
SC = "urban_plume_aging"; PM = ROOT / "runs" / "wp2" / SC / "partmc_s1"; OUT = ROOT / "runs" / "e5_eabs"; OUT.mkdir(parents=True, exist_ok=True)
PROG = {"mixing": 2, "perclass": True, "prognostic_moments": True, "dryshell_selector": 0}
DEF = {"mixing": 2, "perclass": True, "prognostic_moments": False, "dryshell_selector": 0}
SEC = [0.0, 0.5, 1.0]
CONF = {"bc2": {"preset": "custom", "external": True, "primary_fraction_bounds": [0.0, 0.5, 1.0], "secondary_fraction_bounds": SEC},
        "bc4": {"preset": "custom", "external": True, "primary_fraction_bounds": [0.0, 0.25, 0.5, 0.75, 1.0], "secondary_fraction_bounds": SEC},
        "bc7": {"preset": "custom", "external": True, "primary_fraction_bounds": [0.0, 0.1, 0.2, 0.35, 0.5, 0.65, 0.8, 1.0], "secondary_fraction_bounds": SEC},
        "bc10": {"preset": "custom", "external": True, "primary_fraction_bounds": [round(0.1 * i, 2) for i in range(11)], "secondary_fraction_bounds": SEC},
        "NC21": {"preset": "custom", "external": True, "primary_fraction_bounds": [0.0, 0.01, 0.1, 0.25, 0.5, 1.0], "secondary_fraction_bounds": SEC},
        "NC9BC": {"preset": "NC9BC"}}
def ev(d):
    r = load_run(d); k = len(r["results"]["time_s"]) - 1; p = populations(r, k)
    return {"n_class": len(r["grid"]["class_labels"]), "eabs_prog": optics.compute(p, PROG)["bc_abs_enhancement"],
            "eabs_default": optics.compute(p, DEF)["bc_abs_enhancement"], "chi": chi_sections(r["results"]["mass_ugm3"][k], r["results"]["number_m3"][k], r["grid"]["species"]["aerosol"])}
def scram(k):
    c = copy.deepcopy(scn.get(SC)["case"]); c["mixing"] = CONF[k]; c["postprocess"] = {}
    d = OUT / f"scram_{k}"; run_case(c, d); return "scram_" + k, ev(d)
def pmc(k):
    d = OUT / f"pmc_{k}"; project(PM, d, CONF[k]); return "pmc_" + k, ev(d)
res = {}
to = truth_optics(PM, every=48, opts=PROG); res["truth"] = {"time_s": to["time_s"], "eabs": to["bc_abs_enhancement"]}
with ThreadPoolExecutor(12) as ex:
    for fu in [ex.submit(f, k) for k in CONF for f in (scram, pmc)]:
        try:
            k, v = fu.result(); res[k] = v; print(k, v, flush=True)
        except Exception as e:
            print("ERROR", str(e)[-300:], flush=True)
(OUT / "e5.json").write_text(json.dumps(res, indent=1, default=float)); print("=== e5 done")
