"""E2: dust, 8x size sections, six classes: integrator x sub-step scan (plus one-class replicas of r3)."""
import copy, json, sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scrambox import scenarios as scn
from scrambox.case import set_path
from scrambox.post import run_all
from scrambox.run import load_run, run_case
OUT = Path(__file__).resolve().parents[1] / "runs" / "e2_dust_solvers"; OUT.mkdir(parents=True, exist_ok=True)
V = {}
for mix in ("NC6",):
    for solver, subs in (("ros2", (0.5, 1.0, 2.0, 10.0, None)), ("etr", (0.25, 0.5, 1.0, None)), ("euler", (0.5, 1.0, 2.0))):
        for s in subs:
            V[f"dust_x8_{mix}_{solver}_sub{'default' if s is None else f'{s:g}'}"] = (mix, solver, s)
V["dust_x8_NC1_ros2_sub1"] = ("NC1", "ros2", 1.0); V["dust_x8_NC1_etr_sub1"] = ("NC1", "etr", 1.0)
def go(name):
    mix, solver, s = V[name]
    c = copy.deepcopy(scn.get("dust_crustal")["case"])
    c["mixing"] = {"preset": mix}
    c["size_bins"] = {"mode": "log", "n": 14 * 8, "dmin_um": 0.01, "dmax_um": 10.0}
    c["postprocess"] = {"ccn": {"enabled": True, "every": 2}}
    c = set_path(c, "processes.dynamic_solver", solver)
    if s is not None:
        c = set_path(c, "time.substep_s", s)
    d = OUT / name
    try:
        run_case(c, d, timeout_s=40000); r = load_run(d); run_all(r, c["postprocess"], d); r = load_run(d)
        g = {str(x["S_pct"]): round(x["ccn_cm3"], 1) for x in r["post"]["ccn"]["resolved"][-1]["spectrum"]}
        return name, {"N": round(float(r["results"]["number_m3"][-1].sum() / 1e6), 1), "ccn0.1": g.get("0.1"), "ccn0.2": g.get("0.2"),
                      "cpu_s": round(r["summary"].get("cpu_s") or 0, 1), "kernel": r["summary"].get("kernel")}
    except Exception as e:
        return name, {"error": str(e)[-300:]}
res = {}
with ThreadPoolExecutor(int(sys.argv[1]) if len(sys.argv) > 1 else 8) as ex:
    for name, v in ex.map(go, V):
        res[name] = v; print(name, v, flush=True)
        (OUT / "e2.json").write_text(json.dumps(res, indent=1))
print("=== e2 done")
