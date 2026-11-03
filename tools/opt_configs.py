"""Design space for the size x species-group x composition-class optimisation (letter v4).

Three degrees of freedom of SCRAM:
  size grid      : uniform log grids (x1, x2, x4, x8 of the 16-section base) and grids refined
                   only between 0.03 and 0.3 um (CCN critical diameters, nitrate-growth tail)
  species groups : which species define the classification dimensions (case species.groups)
  class bounds   : mass-fraction bounds per group (mixing preset / custom / explicit)
A configuration is CLASS_SCHEME@SIZE_GRID.  mixing_and_overrides(name) -> (mixing, overrides)
for scrambox.partmc.project(run_dir, out_dir, mixing, overrides) or a case copy for run_case.
"""
import math

import numpy as np

DMIN, DMAX, NBASE = 0.005, 10.0, 16
DLOG = math.log10(DMAX / DMIN) / NBASE            # base spacing, 0.206 decades


def _logspace(lo, hi, n):
    return list(np.exp(np.linspace(math.log(lo), math.log(hi), n + 1)))


def refined_bounds(k, lo=0.03, hi=0.3):
    """Base spacing outside [lo, hi], k-fold finer inside."""
    n_lo = max(1, round(math.log10(lo / DMIN) / DLOG))
    n_mid = max(1, round(k * math.log10(hi / lo) / DLOG))
    n_hi = max(1, round(math.log10(DMAX / hi) / DLOG))
    b = _logspace(DMIN, lo, n_lo)[:-1] + _logspace(lo, hi, n_mid)[:-1] + _logspace(hi, DMAX, n_hi)
    return [float(x) for x in b]


SIZE_GRIDS = {
    "U1": {"mode": "log", "n": NBASE, "dmin_um": DMIN, "dmax_um": DMAX},
    "U2": {"mode": "log", "n": 2 * NBASE, "dmin_um": DMIN, "dmax_um": DMAX},
    "U4": {"mode": "log", "n": 4 * NBASE, "dmin_um": DMIN, "dmax_um": DMAX},
    "U8": {"mode": "log", "n": 8 * NBASE, "dmin_um": DMIN, "dmax_um": DMAX},
    "L2": {"mode": "explicit", "bounds_um": refined_bounds(2)},
    "L4": {"mode": "explicit", "bounds_um": refined_bounds(4)},
    "L8": {"mode": "explicit", "bounds_um": refined_bounds(8)},
}

DRY = ["PMD", "PBC", "PNA", "PSO4", "PNH4", "PNO3", "PHCL",
       "PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P", "PPOAlP"]
SOA = ["PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P"]


def _groups(order, members):
    gof = {s: g for g, ss in members.items() for s in ss}
    assert sorted(gof) == sorted(DRY), sorted(set(DRY) ^ set(gof))
    return {"groups": order, "group_of": gof}


GROUPS = {
    # G0: production grouping. SUF = sulfate only; NH4, NO3 and all SOA in the remainder SED.
    "G0": None,
    # G1: same six dimensions, inorganic secondary salts together, SOA as the remainder
    "G1": _groups(["BC", "DU", "SLT", "INO", "POA", "SOA"],
                  {"BC": ["PBC"], "DU": ["PMD"], "SLT": ["PNA", "PHCL"], "INO": ["PSO4", "PNH4", "PNO3"],
                   "POA": ["PPOAlP"], "SOA": SOA}),
    # G2: hygroscopicity only: insoluble/weakly hygroscopic, salts, secondary organics
    "G2": _groups(["INS", "SALT", "ORG"],
                  {"INS": ["PBC", "PMD", "PPOAlP"], "SALT": ["PSO4", "PNH4", "PNO3", "PNA", "PHCL"], "ORG": SOA}),
    # G3: G2 with BC kept apart (absorption)
    "G3": _groups(["BC", "INS", "SALT", "ORG"],
                  {"BC": ["PBC"], "INS": ["PMD", "PPOAlP"], "SALT": ["PSO4", "PNH4", "PNO3", "PNA", "PHCL"], "ORG": SOA}),
    # G4: hygroscopicity groups with sea salt kept apart from the secondary salts
    "G4": _groups(["INS", "SEC", "SS", "ORG"],
                  {"INS": ["PBC", "PMD", "PPOAlP"], "SEC": ["PSO4", "PNH4", "PNO3"], "SS": ["PNA", "PHCL"], "ORG": SOA}),
    # G5: G4 with BC kept apart (absorption)
    "G5": _groups(["BC", "INS", "SEC", "SS", "ORG"],
                  {"BC": ["PBC"], "INS": ["PMD", "PPOAlP"], "SEC": ["PSO4", "PNH4", "PNO3"], "SS": ["PNA", "PHCL"], "ORG": SOA}),
}


def _custom(prim, sec):
    return {"preset": "custom", "external": True, "primary_fraction_bounds": prim, "secondary_fraction_bounds": sec}


B05 = [0.0, 0.5, 1.0]
B3 = [0.0, 0.1, 0.5, 1.0]
B4 = [0.0, 0.2, 0.5, 0.8, 1.0]

# class scheme -> (groups key, mixing)
CLASS_SCHEMES = {
    "NC1":      ("G0", {"preset": "NC1"}),
    "NC6":      ("G0", {"preset": "NC6"}),
    "NC11":     ("G0", {"preset": "NC11"}),
    "NC9BC":    ("G0", {"preset": "NC9BC"}),
    "G1_6":     ("G1", _custom(B05, B05)),
    "G1_9BC":   ("G1", {"preset": "NC9BC"}),
    "G2_3":     ("G2", _custom(B05, B05)),
    "G2_5":     ("G2", _custom(B3, B05)),
    "G2_8":     ("G2", _custom(B3, B3)),
    "G2_10":    ("G2", _custom(B4, B4)),
    "G3_4":     ("G3", _custom(B05, B05)),
    "G3_20":    ("G3", _custom(B3, B3)),
    "G4_4":     ("G4", _custom(B05, B05)),
    "G5_5":     ("G5", _custom(B05, B05)),
}


def mixing_and_overrides(name):
    """'<class scheme>@<size grid>' -> (mixing, overrides for scrambox.case.set_path)."""
    cs, sg = name.split("@")
    gk, mixing = CLASS_SCHEMES[cs]
    ov = {"size_bins": dict(SIZE_GRIDS[sg])}
    if GROUPS[gk] is not None:
        ov["species.groups"] = GROUPS[gk]
    return dict(mixing), ov


def all_configs(schemes=None, grids=None):
    return [f"{c}@{g}" for c in (schemes or CLASS_SCHEMES) for g in (grids or SIZE_GRIDS)]
