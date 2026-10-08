"""Size sections, composition classes and the mapping of log-normal modes onto them.

Composition classes follow CESM-SCRAM exactly (carma_model_mod, CARMA_InitializeModel):
group 1 is the primary group (Primary_Fraction_bounds), groups 2..G-1 are secondary
groups (Secondary_Fraction_bounds), group G is the remainder [0, 1]; a class is one
fraction section per primary/secondary group, kept when the lower bounds sum to < 1.
Section index convention (SCRAM): j = N_frac*(b-1) + f, i.e. bin-major, class stride 1.
"""
from __future__ import annotations

import itertools
import math

import numpy as np

PRESETS = {
    # (external mixing, primary bounds, secondary bounds) as in the CESM-SCRAM profiles
    "NC1": (False, [0.0, 0.5, 1.0], [0.0, 0.5, 1.0]),
    "NC6": (True, [0.0, 0.5, 1.0], [0.0, 0.5, 1.0]),
    "NC11": (True, [0.0, 0.5, 0.9, 1.0], [0.0, 0.5, 0.9, 1.0]),
}


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


def composition_classes(n_groups: int, mixing: dict) -> list[list[tuple[float, float]]]:
    """List of classes; each class is a list of (lower, upper) per group (1..n_groups)."""
    preset = mixing.get("preset", "custom")
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
    """0-based class index of a dry composition (fraction per group, sums to 1)."""
    best, best_d = 0, float("inf")
    for c, cls in enumerate(classes):
        ok = True
        d = 0.0
        for f, (lo, hi) in zip(group_fractions, cls):
            inside = (lo <= f < hi) or (hi == 1.0 and f <= hi + 1e-12)
            if not inside:
                ok = False
                d += min(abs(f - lo), abs(f - hi))
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
