#!/bin/bash
# Post-processing checks: Mie reference values, then optics / CCN / hetchem on the
# mixing-state scenario runs (runs/scn_*_v23 from tests/scenarios_check.sh).
cd /data/shupengz/Model/scram-boxapp2
python3 - <<'EOF'
import math, time
import numpy as np
from scrambox.post.mie import mie_q
# Bohren & Huffman (1983) Appendix A test: m = 1.55, x = 5.213 -> Qext = Qsca = 3.10543
qe, qs, g = mie_q([5.213], [1.55 + 0j])
print(f"BHMIE reference  Qext {qe[0]:.5f} (3.10543)  Qsca {qs[0]:.5f}  g {g[0]:.5f}")
# Rayleigh limit for an absorbing sphere: Qabs -> 4 x Im((m2-1)/(m2+2))
m = 1.85 + 0.71j; x = 1e-3
qe, qs, g = mie_q([x], [m])
print(f"Rayleigh Qabs {qe[0]-qs[0]:.6e}  expected {4*x*((m*m-1)/(m*m+2)).imag:.6e}")
from scrambox.post.common import mam4_wet_radius_m
for rh in (0.3, 0.6, 0.9, 0.98):
    print(f"MAM4 growth factor kappa 0.507, r 50 nm, RH {rh}: {mam4_wet_radius_m(50e-9, 0.507, rh)/50e-9:.3f}")

from scrambox.run import load_run
from scrambox.post import optics, ccn, hetchem
from scrambox.post.common import populations
for sid in ["urban_plume_aging", "seasalt_nitrate", "dust_crustal"]:
    for k in ["v1", "v23"]:
        try:
            r = load_run(f"runs/scn_{sid}_{k}")
        except FileNotFoundError:
            continue
        p = populations(r, -1)
        t0 = time.time()
        o = {m: optics.compute(p, {"mixing": m}) for m in (0, 1, 2, 3)}
        om = optics.compute(p.internal_mean(), {"mixing": 2})
        t1 = time.time()
        c = ccn.compute(p); cm = ccn.compute(p.internal_mean())
        h = hetchem.compute(p)
        t2 = time.time()
        print(f"\n{sid} {k}  t={r['results']['time_s'][-1]/3600:.1f} h  T={p.T:.1f} RH={p.RH:.2f}  "
              f"N={p.number.sum()/1e6:.0f} cm-3  (optics {t1-t0:.1f} s, ccn+het {t2-t1:.1f} s)")
        print("  ext550 Mm-1 by mixing 0/1/2/3: " + "  ".join(f"{o[m]['ext550_Mm']:.1f}" for m in range(4)) +
              f"   internal-mean(2): {om['ext550_Mm']:.1f}")
        print("  ssa550 by mixing 0/1/2/3:      " + "  ".join(f"{o[m]['ssa550']:.3f}" for m in range(4)) +
              f"   internal-mean(2): {om['ssa550']:.3f}")
        print(f"  BC abs enhancement (mixing 1/2): {o[1]['bc_abs_enhancement']}  {o[2]['bc_abs_enhancement']}")
        print("  CCN(S%) resolved: " + "  ".join(f"{s['S_pct']}:{s['ccn_cm3']:.0f}" for s in c["spectrum"]))
        print("  CCN(S%) int-mean: " + "  ".join(f"{s['S_pct']}:{s['ccn_cm3']:.0f}" for s in cm["spectrum"]))
        print("  Nd(w) resolved:   " + "  ".join(f"{a['w_m_s']}:{a['nd_cm3']:.0f}(s={a['smax_pct'] or 0:.3f}%)" for a in c["activation"]))
        print("  Nd(w) int-mean:   " + "  ".join(f"{a['w_m_s']}:{a['nd_cm3']:.0f}(s={a['smax_pct'] or 0:.3f}%)" for a in cm["activation"]))
        print(f"  SAD {h['sad_um2_cm3']:.0f} um2/cm3")
        for gas in ("N2O5", "HO2"):
            for lab in ("resolved", "internal_mean"):
                rr = h["rules"][gas][lab]
                print(f"  {gas:4s} {lab:13s} " + "  ".join(f"{rule}: k={v['k_s1']:.2e} g_eff={v['gamma_eff']:.3g}" for rule, v in rr.items()))
EOF
