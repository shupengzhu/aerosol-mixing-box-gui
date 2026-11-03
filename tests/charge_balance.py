"""Ion balance per size section of a finished run (dust + crustal ions check).

Cations: NH4+, Na+ (PNA, plus dust Na), Ca2+, K+, Mg2+ diagnosed from PMD and PNA with the
CMAQ fractions the v2.3 kernel uses (kernel/data/speciation_cmaq_aero.json).
Anions: SO4 2-, NO3-, Cl-. Reported in neq m-3 and as the cation/anion ratio.
Usage: python3 tests/charge_balance.py RUN_DIR [RUN_DIR ...]
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scrambox.run import load_run  # noqa: E402

SPEC = json.loads((Path(__file__).resolve().parent.parent / "kernel/data/speciation_cmaq_aero.json").read_text())
MW = {"PNH4": 18.04, "PNA": 22.99, "PSO4": 96.06, "PNO3": 62.0, "PHCL": 35.45, "Ca": 40.08, "K": 39.10, "Mg": 24.31, "Na": 22.99}
CHG = {"PNH4": 1, "PNA": 1, "PSO4": 2, "PNO3": 1, "PHCL": 1, "Ca": 2, "K": 1, "Mg": 2, "Na": 1}


def balance(run_dir, it=-1, crustal=True):
    r = load_run(run_dir)
    sp, b = r["grid"]["species"]["aerosol"], np.array(r["grid"]["bounds_um"])
    m = r["results"]["mass_ugm3"][it].sum(axis=1)               # (bin, species), classes summed
    mid = np.sqrt(b[:-1] * b[1:])
    rows = []
    for k in range(m.shape[0]):
        get = lambda s: m[k, sp.index(s)] if s in sp else 0.0
        eq = lambda mass, ion: 1e3 * mass / MW[ion] * CHG[ion]  # ug/m3 -> neq/m3
        cat = eq(get("PNH4"), "PNH4") + eq(get("PNA"), "PNA")
        if crustal:
            cut = "fine" if b[k + 1] <= SPEC["fine_coarse_cut_diameter_um"] + 1e-9 else "coarse"
            fr = SPEC["dust_fraction_of_parent"][cut]
            for ion in ("Ca", "K", "Mg", "Na"):
                cat += eq(get("PMD") * fr[ion], ion)
            # sea-salt Ca/K/Mg re-speciate the PNA sodium equivalents (kernel
            # commit after f86234e), so they add no equivalents; before that fix they were
            # added on top: pass --ss-additive to reproduce the old accounting
            if "--ss-additive" in sys.argv:
                for ion, rto in SPEC["sea_salt_ratio_to_Na"].items():
                    cat += eq(get("PNA") * rto, ion)
        an = eq(get("PSO4"), "PSO4") + eq(get("PNO3"), "PNO3") + eq(get("PHCL"), "PHCL")
        rows.append((mid[k], cat, an, get("PNO3"), get("PH2O")))
    return r, rows


for d in [a for a in sys.argv[1:] if not a.startswith("--")]:
    crustal = "v23" in d
    r, rows = balance(d, crustal=crustal)
    print(f"\n{d}  (crustal cations {'included' if crustal else 'ignored (v1)'}), t = {r['results']['time_s'][-1]/3600:.1f} h")
    print("  D_um      cat_neq   an_neq   cat/an   NO3_ug   H2O_ug")
    for mid, c, a, no3, h2o in rows:
        if c + a > 1e-3:
            print(f"  {mid:7.3f} {c:9.2f} {a:8.2f} {c/a if a > 0 else float('nan'):7.2f} {no3:8.3f} {h2o:8.2f}")
    c = sum(x[1] for x in rows); a = sum(x[2] for x in rows)
    print(f"  total   {c:9.2f} {a:8.2f} {c/a:7.2f}")
