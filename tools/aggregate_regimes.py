"""Aggregate the regime ensemble into one compact JSON (end-of-run and time-mean errors, MAM4,
numerics, E_abs with prognostic-moment optics, droplet-number chain). Usage: python tools/aggregate_regimes.py out.json"""
import json, sys, glob
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from scrambox.post import optics
from scrambox.post.common import populations
from scrambox.run import load_run
OPTS = {"mixing": 2, "perclass": True, "prognostic_moments": True, "dryshell_selector": 0}
SUB = sys.argv[2] if len(sys.argv) > 2 else "regime"
IDX = sys.argv[3] if len(sys.argv) > 3 else str(Path.home() / "grl_box" / "regime_scen" / "regime_index.json")
idx = {r["id"]: r for r in json.load(open(IDX))}
def eabs(d):
    try:
        r = load_run(d); o = optics.compute(populations(r, len(r["results"]["time_s"]) - 1), OPTS)
        return o["bc_abs_enhancement"], o["abs550_Mm"]
    except Exception as e:
        return None, None
out = []
for S in sorted(p for p in (ROOT / "runs" / SUB).iterdir() if p.is_dir() and p.name in idx):
    sid = S.name; row = {"id": sid, **{k: v for k, v in idx[sid].items() if k not in ("gas_src",)}}
    try:
        sm = json.loads((S / "summary.json").read_text())
    except Exception:
        row["error"] = "no summary"; out.append(row); continue
    T = sm["truth"]["s1"]; tt = np.array(T["time_s"]) / 3600
    row["truth"] = {"ccn": {s: T["ccn_cm3"][s][-1] for s in ("0.1", "0.2", "0.5")}, "chi": T["chi"][-1], "N": T["n_total_cm3"][-1]}
    for nc in ("NC1", "NC6", "NC11"):
        for lab, key in (("rep", f"partmc_s1_{nc}"), ("sim", f"scram_{nc}")):
            R = sm["runs"][key]; th = np.array(R["post_time_h"])
            e = {"ccn": {s: R["ccn"][s][-1] for s in ("0.1", "0.2", "0.5")}, "chi": R["chi_sections"][-1], "PM": R["PM_dry"][-1],
                 "bc": R["species"].get("PBC", [0])[-1]}
            tm = {}
            for s in ("0.1", "0.2", "0.5"):
                tru = np.interp(th, tt, np.array(T["ccn_cm3"][s])); m = (th >= 0.999) & (tru > 1)
                rel = np.array(R["ccn"][s])[m] / tru[m] - 1
                tm[s] = {"mean_abs": float(np.mean(np.abs(rel))) if m.any() else None, "max_abs": float(np.max(np.abs(rel))) if m.any() else None}
            e["time"] = tm
            kb = R.get("k_n2o5_bt"); kc = R.get("k_n2o5_cesm")
            e["k_bt"] = kb[-1] if kb else None; e["k_cesm"] = kc[-1] if kc else None
            row[f"{lab}_{nc}"] = e
    for lab, d in (("sim_NC6", S / "scram_NC6"), ("rep_NC6", S / "partmc_s1_NC6"), ("sim_NC1", S / "scram_NC1"), ("rep_NC1", S / "partmc_s1_NC1")):
        row[lab]["eabs"], row[lab]["abs550"] = eabs(d)
    try:
        h = json.loads((S / "het_truth.json").read_text()); row["het"] = h
    except Exception:
        row["het"] = None
    try:
        o = json.loads((S / "optics.json").read_text()); row["truth"]["eabs"] = o["truth"]["1"]["bc_abs_enhancement"][-1] if "1" in o["truth"] else o["truth"][list(o["truth"])[0]]["bc_abs_enhancement"][-1]
    except Exception:
        row["truth"]["eabs"] = None
    try:
        st = json.loads((S / "stage2.json").read_text())
        row["mam4"] = {"ccn": st["mam4"]["ccn"]}
        row["activation"] = st["activation"]
        num = st.get("numerics", {})
        for lab in num:
            for who in ("scram", "pmc"):
                d = S / "stage2" / f"{who}_{lab}"; ea, ab = eabs(d)
                num[lab][who]["eabs_prog"] = ea
        row["numerics"] = num
    except Exception as e:
        row["stage2_error"] = str(e)[-200:]
    out.append(row)
json.dump(out, open(sys.argv[1], "w"), default=float)
print(len(out), "regimes;", sum(1 for r in out if "stage2_error" not in r and "error" not in r), "complete")
