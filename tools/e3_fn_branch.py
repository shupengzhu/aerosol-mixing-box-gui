"""E3: FN population-splitting branch (discriminant D) on the 732 synthetic cases, scrambox constants."""
import json, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent)); sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parcel_activation import coefficients
A = json.load(open(sys.argv[1]))
T, P = 283.15, 9.0e4
alpha, gamma, G, Ak = coefficients(T, P)
rows = []
for r in A:
    if r["truth"] <= 1: continue
    D = r["smax_FN"] ** 4 - 16 * Ak ** 2 * alpha * r["w"] / (9 * G)
    rows.append((D >= 0, r["FN"] / r["truth"], r["current"] / r["truth"], r["count_at_true_smax"] / r["truth"], r["smax_FN"] / r["smax_true"], r["smax_current"] / r["smax_true"], r["w"], r["label"]))
def st(x):
    x = np.asarray(x, float); return {"median": round(float(np.median(x)), 3), "p10": round(float(np.percentile(x, 10)), 3), "p90": round(float(np.percentile(x, 90)), 3), "rms_log10": round(float(np.sqrt(np.mean(np.log10(x) ** 2))), 3)}
out = {"constants": {"T": T, "P": P, "alpha": alpha, "G": G, "A": Ak}}
for flag, lab in ((True, "D>=0 analytic"), (False, "D<0 empirical")):
    sel = [r for r in rows if r[0] == flag]
    out[lab] = {"n": len(sel), "FN": st([r[1] for r in sel]), "ARG": st([r[2] for r in sel]), "count_exact_smax": st([r[3] for r in sel]),
                "smax_FN_ratio": st([r[4] for r in sel]), "smax_ARG_ratio": st([r[5] for r in sel]),
                "by_w": {str(w): sum(1 for r in sel if r[6] == w) for w in sorted({r[6] for r in rows})},
                "by_type": {t: sum(1 for r in sel if r[7].startswith(t)) for t in ("uni", "bi", "ext")}}
out["FN_over_count_all"] = st([r[1] / r[3] for r in rows])
json.dump(out, open(sys.argv[2], "w"), indent=1); print(json.dumps(out, indent=1))
