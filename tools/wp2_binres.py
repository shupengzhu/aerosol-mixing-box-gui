"""Is the dust CCN(0.1%) gap SCRAM vs PartMC sectional numerical diffusion? SCRAM NC1 with the
scenario's bins refined x2 / x4 (same range) and alternative redistribution schemes."""
import copy, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import numpy as np
sys.path.insert(0, "/data/shupengz/Model/scram-boxapp2")
sys.path.insert(0, "/data/shupengz/Model/scram-boxapp2/tools")
from scrambox import scenarios as scn
from scrambox.post import run_all
from scrambox.run import load_run, run_case
from partmc_experiment import series

sc = sys.argv[1] if len(sys.argv) > 1 else "dust_crustal"
OUT = Path("/data/shupengz/Model/scram-boxapp2/runs/wp2_binres") / sc
OUT.mkdir(parents=True, exist_ok=True)
base = copy.deepcopy(scn.get(sc)["case"])
base["mixing"] = {"preset": "NC1"}
base["postprocess"] = {"ccn": {"enabled": True, "every": 6}}
b0 = np.asarray(load_run(Path("/data/shupengz/Model/scram-boxapp2/runs/wp2") / sc / "scram_NC1")["grid"]["bounds_um"])
print("base bins", len(b0) - 1, b0[0], b0[-1])

def refine(f):
    lb = np.log(b0)
    return np.exp(np.concatenate([np.linspace(lb[i], lb[i + 1], f + 1)[:-1] for i in range(len(b0) - 1)] + [[lb[-1]]])).tolist()

variants = {"x1": {"bins": refine(1)}, "x2": {"bins": refine(2)}, "x4": {"bins": refine(4)},
            "x8": {"bins": refine(8)},
            "x1_ppm_dual": {"bins": refine(1), "redist": "core_ppm_dualmoment"},
            "x1_legacy": {"bins": refine(1), "redist": "legacy"}}

if "--refine-only" in sys.argv:
    variants = {k: v for k, v in variants.items() if "redist" not in v}


def go(name):
    v = variants[name]; c = copy.deepcopy(base)
    if "bins" in v:
        c["size_bins"] = {"mode": "explicit", "bounds_um": v["bins"]}
    if "redist" in v:
        c.setdefault("processes", {})["redistribution"] = v["redist"]
    d = OUT / name
    try:
        run_case(c, d)
        r = load_run(d)
        run_all(r, c["postprocess"], d)
        s = series(load_run(d))
        return name, {k: s["ccn"][k][-1] for k in s["ccn"]}, s["N_cm3"][-1], s["PM_dry"][-1], s["post_time_h"][-1]
    except Exception as e:
        return name, f"ERROR {e}", None, None, None

with ThreadPoolExecutor(6) as ex:
    for name, ccn, n, pm, th in ex.map(go, variants):
        print(name, "t=", th, "N", n, "PM", pm, "CCN", ccn, flush=True)
