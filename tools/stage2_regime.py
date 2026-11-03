"""Stage 2 for one regime run directory (runs/regime/<id>): MAM4 static projection, droplet-number
error chain (parcel / ARG / FN), size refinement (NC6 x8) and BC-resolved classes (NC9BC x1).
Usage: python tools/stage2_regime.py <id> [--skip-numerics]   -> runs/regime/<id>/stage2.json"""
import copy, json, math, sys
from pathlib import Path
import numpy as np
T_ = Path(__file__).resolve().parent; ROOT = T_.parent
sys.path.insert(0, str(T_)); sys.path.insert(0, str(ROOT))
from parcel_activation import parcel_nd
from wp2_activation_parcel import last_particles, truth_parcel, mam4_parcel, scram_pops
from wp2_mam4_projection import assign_modes, mam4_metrics
import activation_validation as av
from scrambox import scenarios as scn
from scrambox.partmc import chi_sections, project
from scrambox.post import run_all
from scrambox.post import ccn as C
from scrambox.post.common import populations
from scrambox.run import load_run, run_case
import os
sid = sys.argv[1]; S = ROOT / "runs" / os.environ.get("STAGE2_RUNS", "regime") / sid; W = (0.2, 0.5, 1.0)
out = {"id": sid}
rec = last_particles(S / "partmc_s1")
t, T, P, names, vol, n, kap = rec
mm = mam4_metrics(T, P, names, vol, n, kap, assign_modes(names, vol, kap), True)
out["mam4"] = {"ccn": mm["ccn"], "nd_arg": mm["nd"], "modes": mm["modes"]}
def fn_arg(run_dir, w):
    pop = populations(load_run(run_dir), -1); md = C.modes(pop)
    av.T, av.P = pop.T, pop.P; zeta, pref, A = av.arg_consts(w)
    sl = [{"N": m["n_m3"], "sc": m["sc"], "dlo": m["dlo"], "dhi": m["dhi"], "tilt": m["tilt"], "kappa": m["kappa"], "dmid": m["dmid"]} for m in md]
    s_arg = av.maxsat([s["N"] for s in sl], [s["sc"] for s in sl], 0.5, 1.0, zeta, pref)
    Nf, scf = av.fn_spectrum(sl); st = []; s_fn = av.smax_fn(Nf, scf, w, pop.T, pop.P, stats=st)
    nd = lambda s: sum(C._activated(m, s) for m in md) / 1e6
    return {"arg": nd(s_arg), "arg_smax": s_arg * 100, "fn": nd(s_fn), "fn_smax": s_fn * 100, "fn_D_pos": bool(st[0])}
def parcel_on(run_dir, w):
    pop, nn, rd, kk = scram_pops(run_dir); sm, nd = parcel_nd(pop.T, pop.P, w, nn, rd, kk); return {"parcel": nd / 1e6, "parcel_smax": sm * 100}
act = {}
for w in W:
    sm, nd = truth_parcel(rec, w); row = {"truth": {"nd": nd / 1e6, "smax": sm * 100}}
    sm2, nd2 = mam4_parcel(rec, w); row["mam4_parcel"] = {"nd": nd2 / 1e6, "smax": sm2 * 100}
    for lab, d in (("p_nc6", S / "partmc_s1_NC6"), ("scram_nc6", S / "scram_NC6")):
        row[lab] = {**parcel_on(d, w), **fn_arg(d, w)}
    act[str(w)] = row; print(sid, "w", w, {k: (v.get("nd") or v.get("parcel")) for k, v in row.items()}, flush=True)
out["activation"] = act
if "--skip-numerics" not in sys.argv:
    POST = {"optics": {"enabled": True, "every": 24}, "ccn": {"enabled": True, "every": 24}}
    base = copy.deepcopy(scn.get(sid)["case"])
    def summ(d):
        r = load_run(d); m, nn = r["results"]["mass_ugm3"][-1], r["results"]["number_m3"][-1]; sp = r["grid"]["species"]["aerosol"]
        o = r["post"]["optics"]["resolved"][-1]; cc = {str(x["S_pct"]): x["ccn_cm3"] for x in r["post"]["ccn"]["resolved"][-1]["spectrum"]}
        ib = sp.index("PBC")
        return {"chi": chi_sections(m, nn, sp), "eabs": o["bc_abs_enhancement"], "bc_ugm3": float(m[..., ib].sum()),
                "ccn": {k: cc[k] for k in ("0.1", "0.2", "0.5") if k in cc}, "cpu_s": (r.get("summary") or {}).get("cpu_s")}
    num = {}
    for lab, mix, nb in (("NC6_x8", "NC6", 128), ("NC9BC_x1", "NC9BC", 16), ("NC6_x1", "NC6", 16)):
        bins = {"mode": "log", "n": nb, "dmin_um": 0.005, "dmax_um": 10.0}
        c = copy.deepcopy(base); c["mixing"] = {"preset": mix}; c["size_bins"] = bins; c["postprocess"] = POST
        d = S / "stage2" / f"scram_{lab}"; run_case(c, d, timeout_s=40000); run_all(load_run(d), POST, d)
        dp = S / "stage2" / f"pmc_{lab}"; project(S / "partmc_s1", dp, {"preset": mix}, {"size_bins": bins}); run_all(load_run(dp), POST, dp)
        num[lab] = {"scram": summ(d), "pmc": summ(dp)}; print(sid, lab, num[lab], flush=True)
    out["numerics"] = num
tr = json.loads((S / "partmc_s1" / "truth.json").read_text())
out["truth"] = {"chi": tr["chi"][-1], "ccn": {s: tr["ccn_cm3"][s][-1] for s in ("0.1", "0.2", "0.5")}, "n_total": tr["n_total_cm3"][-1],
                "bc_coating_ratio": tr.get("bc_coating_ratio", [None])[-1]}
(S / "stage2.json").write_text(json.dumps(out, indent=1, default=float)); print("=== stage2 done", sid)
