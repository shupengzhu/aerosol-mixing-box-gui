#!/usr/bin/env python3
"""Convergence of the parcel reference on the urban MAM4 modes (end of run, seed 1):
start supersaturation, tolerance and section count must not change S_max / Nd."""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp  # noqa: F401  (parcel_activation uses it)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import parcel_activation as pa                                   # noqa: E402
from wp2_activation_parcel import last_particles, lognormal_sections   # noqa: E402
from wp2_mam4_projection import assign_modes, mam4_metrics      # noqa: E402

root = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/wp2")
rec = last_particles(root / "urban_plume_aging" / "partmc_s1")
t, T, P, names, vol, n, kap = rec
modes = mam4_metrics(T, P, names, vol, n, kap, assign_modes(names, vol, kap), True)["modes"]
print("T", T, "P", P, [(m["mode"], round(m["N_cm3"]), round(m["dgn_um"], 4), round(m["kappa"], 3)) for m in modes])
for nsec in (40, 120, 400):
    ns, rs, ks = [], [], []
    for m in modes:
        a, r = lognormal_sections(m["N_cm3"] * 1e6, m["dgn_um"], m["sigma"], nsec)
        ns.append(a); rs.append(r); ks.append(np.full(len(a), m["kappa"]))
    N, R, K = map(np.concatenate, (ns, rs, ks))
    for S0 in (-0.005, -0.02):
        for w in (0.2, 1.0):
            sm, nd = pa.parcel_nd(T, P, w, N, R, K, S0=S0)
            print(f"nsec {nsec:4d} S0 {S0:+.3f} w {w}: smax {sm * 100:.4f} %  Nd {nd / 1e6:.1f}")
# tighter tolerance via monkeypatch
orig = pa.solve_ivp
def tight(*a, **k):
    k["rtol"] = 1e-9
    k["atol"] = k["atol"] * 1e-3
    k["max_step"] = k["max_step"] / 10
    return orig(*a, **k)
pa.solve_ivp = tight
for w in (0.2, 1.0):
    sm, nd = pa.parcel_nd(T, P, w, N, R, K)
    print(f"tight tol nsec 400 w {w}: smax {sm * 100:.4f} %  Nd {nd / 1e6:.1f}")
