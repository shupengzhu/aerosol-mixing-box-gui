"""Size sections, composition classes and the mapping of log-normal modes onto them.

Composition classes follow SCRAM exactly (carma_model_mod, CARMA_InitializeModel):
group 1 is the primary group (Primary_Fraction_bounds), groups 2..G-1 are secondary
groups (Secondary_Fraction_bounds), group G is the remainder [0, 1]; a class is one
fraction section per primary/secondary group, kept when the lower bounds sum to < 1.
Section index convention (SCRAM): j = N_frac*(b-1) + f, i.e. bin-major, class stride 1.
"""
from __future__ import annotations

import itertools
import math
import os

import numpy as np

PRESETS = {
    # (external mixing, primary bounds, secondary bounds) as in the SCRAM profiles
    "NC1": (False, [0.0, 0.5, 1.0], [0.0, 0.5, 1.0]),
    "NC6": (True, [0.0, 0.5, 1.0], [0.0, 0.5, 1.0]),
    "NC11": (True, [0.0, 0.5, 0.9, 1.0], [0.0, 0.5, 0.9, 1.0]),
    # BC resolved below 0.5 (WP2 2026-10-09: E_abs within ~1 % of PartMC; ~9x the cost of NC6)
    "NC21": (True, [0.0, 0.01, 0.1, 0.25, 0.5, 1.0], [0.0, 0.5, 1.0]),
}

# Presets that are an explicit class list rather than a product of fraction bounds.
# NC9BC: BC-conditional refinement. BC-free particles (BC <= 0.01) keep the NC6 split of the
# other groups; BC-containing particles are classified by their BC fraction only
# (0.01-0.1, 0.1-0.25, 0.25-0.5, 0.5-1), whatever their coating is made of.
BC_CONDITIONAL = {"NC9BC": ([0.0, 0.01, 0.1, 0.25, 0.5, 1.0], [0.0, 0.5, 1.0])}


def _bc_conditional(n_groups: int, bc: list[float], sec: list[float]) -> list[list[tuple[float, float]]]:
    free = composition_classes(n_groups, {"preset": "custom", "external": True,
                                          "primary_fraction_bounds": [bc[0], bc[1], 1.0],
                                          "secondary_fraction_bounds": sec})
    classes = [c for c in free if c[0] == (bc[0], bc[1])]
    for lo, hi in zip(bc[1:-1], bc[2:]):
        classes.append([(lo, hi)] + [(0.0, 1.0)] * (n_groups - 1))
    return classes


def size_bounds_um(spec: dict) -> list[float]:
    if spec.get("mode", "log") == "explicit":
        b = [float(x) for x in spec["bounds_um"]]
    else:
        n, lo, hi = int(spec["n"]), float(spec["dmin_um"]), float(spec["dmax_um"])
        if n < 1 or lo <= 0 or hi <= lo:
            raise ValueError("size_bins: need n >= 1 and 0 < dmin_um < dmax_um")
        b = list(np.exp(np.linspace(math.log(lo), math.log(hi), n + 1)))
    if len(b) < 2 or any(b[i + 1] <= b[i] for i in range(len(b) - 1)) or b[0] <= 0:
        raise ValueError("size_bins: bounds must be positive and strictly increasing")
    return b


def explicit_classes(n_groups: int, rows) -> list[list[tuple[float, float]]]:
    """Classes given one by one (mixing preset "explicit", cfg block Class_bounds): each row
    has (lo, hi) for the first n_groups-1 groups, or for all groups with the last [0, 1]."""
    classes = []
    for r, row in enumerate(rows, 1):
        pairs = [tuple(float(x) for x in p) for p in row]
        if any(len(p) != 2 for p in pairs) or len(pairs) not in (n_groups - 1, n_groups):
            raise ValueError(f"class {r}: needs [lo, hi] for each of the first {n_groups - 1} groups")
        if len(pairs) == n_groups - 1:
            pairs.append((0.0, 1.0))
        classes.append(pairs)
    if not classes:
        raise ValueError("explicit classes: the list is empty")
    return classes


def check_partition(classes, groups: list[str] | None = None) -> list[str]:
    """Closure test of a class list, identical to the kernel's scram_check_class_partition.
    Membership per group: [0, hi] when lo = 0, (lo, hi] otherwise; dimensions are the groups
    except the last one (the remainder, which must be [0, 1]). The breakpoints of all
    intervals cut every dimension into elementary intervals with the same rule, so an
    elementary cell is either inside or outside each class; it holds compositions iff the sum
    of its lower ends is < 1. Every such cell must lie in exactly one class (exact, no
    sampling). Raises ValueError naming the first region in no class or in several classes;
    returns notes on classes that contain no composition."""
    tol = 1e-12
    if not classes:
        raise ValueError("composition classes: the list is empty")
    ng = len(classes[0])
    nd = ng - 1
    names = groups or [f"group {g + 1}" for g in range(ng)]
    for c, cls in enumerate(classes, 1):
        if len(cls) != ng:
            raise ValueError(f"class {c}: wrong number of groups")
        for g, (lo, hi) in enumerate(cls):
            if not (0.0 <= lo < hi <= 1.0) or hi - lo <= tol:
                raise ValueError(f"class {c}, {names[g]}: need 0 <= lo < hi <= 1")
        if abs(cls[-1][0]) > tol or abs(cls[-1][1] - 1.0) > tol:
            raise ValueError(f"class {c}: the last group ({names[-1]}, remainder) must be [0, 1]")
    bps = []
    for d in range(nd):
        v = sorted({0.0, 1.0} | {x for cls in classes for x in cls[d]})
        u = [v[0]]
        for x in v[1:]:
            if x - u[-1] > tol:
                u.append(x)
        bps.append(u)
    used = [False] * len(classes)
    for idx in itertools.product(*[range(len(b) - 1) for b in bps]):
        # A cell intersects the composition simplex whenever its lower-corner
        # sum is below one by the existing tolerance. Its midpoint may be outside while a triangular
        # part of the cell is feasible, so midpoint sampling misses real gaps.
        if sum(bps[d][i] for d, i in enumerate(idx)) >= 1.0 - tol:
            continue
        hits = [c for c, cls in enumerate(classes)
                if all(cls[d][0] <= bps[d][i] + tol and bps[d][i + 1] <= cls[d][1] + tol for d, i in enumerate(idx))]
        if len(hits) != 1:
            region = ", ".join(f"{names[d]} {'[' if i == 0 else '('}{bps[d][i]:g}, {bps[d][i + 1]:g}]"
                               for d, i in enumerate(idx) if (bps[d][i], bps[d][i + 1]) != (0.0, 1.0))
            what = "belong to no class" if not hits else f"belong to {len(hits)} classes ({', '.join(f'C{c + 1}' for c in hits)})"
            raise ValueError(f"composition classes are not closed: compositions with {region or 'any composition'} {what}")
        used[hits[0]] = True
    return [f"class C{c + 1} contains no composition (its lower bounds sum to >= 1)" for c, u in enumerate(used) if not u]


def composition_classes(n_groups: int, mixing: dict, groups: list[str] | None = None) -> list[list[tuple[float, float]]]:
    """List of classes; each class is a list of (lower, upper) per group (1..n_groups)."""
    preset = mixing.get("preset", "custom")
    if preset == "explicit" and n_groups > 1:
        classes = explicit_classes(n_groups, mixing.get("classes") or [])
        if not os.environ.get("SCRAMBOX_SKIP_CLASS_CHECK"):      # kernel-check tests only
            check_partition(classes, groups)
        return classes
    if preset in BC_CONDITIONAL and n_groups > 2:
        classes = _bc_conditional(n_groups, *BC_CONDITIONAL[preset])
        check_partition(classes, groups)
        return classes
    if preset in PRESETS:
        external, prim, sec = PRESETS[preset]
    else:
        external = bool(mixing.get("external", True))
        prim = [float(x) for x in mixing.get("primary_fraction_bounds", [0.0, 1.0])]
        sec = [float(x) for x in mixing.get("secondary_fraction_bounds", [0.0, 1.0])]
    for name, b in (("primary", prim), ("secondary", sec)):
        if b[0] != 0.0 or b[-1] != 1.0 or any(b[i + 1] <= b[i] for i in range(len(b) - 1)):
            raise ValueError(f"{name}_fraction_bounds must increase from 0 to 1")
    if not external or n_groups <= 1:
        return [[(0.0, 1.0)] * n_groups]
    classes = []
    if n_groups == 2:
        for i in range(len(prim) - 1):
            classes.append([(prim[i], prim[i + 1]), (0.0, 1.0)])
        return classes
    nsec = len(sec) - 1
    for i in range(len(prim) - 1):
        # SCRAM's counter runs with the first secondary group varying fastest
        for combo in itertools.product(range(nsec), repeat=n_groups - 2):
            idx = combo[::-1]
            if prim[i] + sum(sec[k] for k in idx) < 1.0:
                cls = [(prim[i], prim[i + 1])]
                cls += [(sec[k], sec[k + 1]) for k in idx]
                cls.append((0.0, 1.0))
                classes.append(cls)
    if not classes:
        raise ValueError("no valid composition class (check fraction bounds)")
    check_partition(classes, groups)
    return classes


def comp_bounds_array(classes) -> list[float]:
    out = []
    for cls in classes:
        for lo, hi in cls:
            out += [lo, hi]
    return out


def class_label(cls, groups: list[str]) -> str:
    parts = []
    for (lo, hi), g in zip(cls[:-1], groups[:-1]):
        if (lo, hi) != (0.0, 1.0) and lo > 0.0:
            parts.append(f"{g} {lo:g}-{hi:g}")
    if not parts:
        dom = [f"{g}<{hi:g}" for (lo, hi), g in zip(cls[:-1], groups[:-1]) if hi < 1.0]
        return "mixed (" + ", ".join(dom) + ")" if dom else "internal mixture"
    return ", ".join(parts)


def find_class(group_fractions: list[float], classes) -> int:
    """0-based class index of a dry composition (fraction per group, sums to 1), with the
    kernel's rule ([0, hi] when lo = 0, (lo, hi] otherwise; groups except the remainder).
    If no class contains it (rounding at a bound), the nearest class, as the kernel's
    scram_nearest_class."""
    best, best_d = 0, float("inf")
    for c, cls in enumerate(classes):
        ok = True
        d = 0.0
        for f, (lo, hi) in zip(group_fractions[:-1], cls[:-1]):
            if not ((lo <= f <= hi) if lo == 0.0 else (lo < f <= hi)):
                ok = False
            d += max(0.0, lo - f) + max(0.0, f - hi)
        if ok:
            return c
        if d < best_d:
            best, best_d = c, d
    return best


def _phi(x):
    return 0.5 * (1.0 + math.erf(x / math.sqrt(2.0)))


def lognormal_to_bins(number_m3: float, dg_um: float, sigma_g: float, bounds_um: list[float]):
    """Number (#/m3) and dry volume (um3/m3) per section for a log-normal mode;
    also the fractions of number and volume that fall outside the section range."""
    ls = math.log(sigma_g)
    if sigma_g <= 1.0 or dg_um <= 0:
        raise ValueError("log-normal mode needs dg_um > 0 and sigma_g > 1")
    vol_total = number_m3 * math.pi / 6.0 * dg_um ** 3 * math.exp(4.5 * ls * ls)
    nb, vb = [], []
    for lo, hi in zip(bounds_um[:-1], bounds_um[1:]):
        zl, zh = math.log(lo / dg_um) / ls, math.log(hi / dg_um) / ls
        nb.append(number_m3 * (_phi(zh) - _phi(zl)))
        vb.append(vol_total * (_phi(zh - 3.0 * ls) - _phi(zl - 3.0 * ls)))
    out_n = 1.0 - sum(nb) / number_m3 if number_m3 > 0 else 0.0
    out_v = 1.0 - sum(vb) / vol_total if vol_total > 0 else 0.0
    return np.array(nb), np.array(vb), out_n, out_v
