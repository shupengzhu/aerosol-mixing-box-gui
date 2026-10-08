"""Cloud activation of a box state, identical to CESM-SCRAM scram_cloud_aerosol_export.

Every (bin, class) population is one monodisperse mode at the bin reference dry radius
with volume-weighted kappa (PH2O excluded). Koehler critical supersaturation
sc = 2 A sqrt(A / (27 kappa r^3)), A = 2 Mw sigma / (R T rho_w). Sub-bin activation uses
the linear dry-mass tilt reconstructed from the bin's number and dry volume.

  ccn_spectrum   CCN at prescribed supersaturations (the CESM CCN3 path, S = 0.1 %)
  activation     Abdul-Razzak & Ghan (2000) maximum supersaturation for an updraft w
                 (f1 = 0.5, f2 = 1 for bins, ndrop_bam thermodynamics, w >= 0.2 m/s)
"""
from __future__ import annotations

import math

import numpy as np

from .common import Populations, kappa_of

# CAM physconst
R_UNIV, MWH2O, MWDRY = 8.31446261815324e3, 18.016, 28.966
RAIR, RH2O = R_UNIV / MWDRY, R_UNIV / MWH2O
LATVAP, CPAIR, GRAVIT, RHOH2O, TMELT = 2.501e6, 1.00464e3, 9.80616, 1.0e3, 273.15
SURFTEN, PREF = 0.076, 1013.25e2
F1_BIN, F2_BIN = 0.5, 1.0
WSUB_MIN = 0.20
NUMBER_FLOOR_CM3 = 1.0e-12
EPSILO = MWH2O / MWDRY

DEFAULTS = {"supersaturations_pct": [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0],
            "updrafts_m_s": [0.2, 0.5, 1.0, 2.0]}


# ---------------------------------------------------------------- linear mass tilt
def _inv_d3_integral(dlo, dhi, ratio, dcut):
    if dlo <= 0 or dhi <= dlo or dcut >= dhi:
        return 0.0
    ratio = max(ratio, 1e-6)
    span = math.log(dhi / dlo)
    x0 = min(1.0, max(0.0, math.log(dcut / dlo) / span)) if dcut > dlo else 0.0
    a = 3.0 * span
    eh, ec = math.exp(-a), math.exp(-a * x0)
    return span / dlo ** 3 * (ratio * (ec - eh) / a + (1.0 - ratio) * (ec * (a * x0 + 1.0) - eh * (a + 1.0)) / (a * a))


def _number_factor(dlo, dhi, ratio):
    ratio = max(ratio, 1e-6)
    span = math.log(dhi / dlo)
    jvol = span * 0.5 * (1.0 + ratio)
    jnum = _inv_d3_integral(dlo, dhi, ratio, dlo)
    return (6.0 / math.pi) * jnum / jvol if jvol > 0 and jnum > 0 else 0.0


def tilt_ratio(dlo, dhi, total_num, total_dryvol):
    """scram_linear_mass_tilt_ratio_from_num_volume (bisection in log ratio)."""
    if dlo <= 0 or dhi <= dlo or total_num <= 0 or total_dryvol <= 0:
        return 1.0
    target = total_num / total_dryvol
    lo = 1e-6
    f_lo = _number_factor(dlo, dhi, lo)
    if f_lo <= 0:
        return 1.0
    if target <= f_lo:
        return lo
    hi = 1.0
    f_hi = _number_factor(dlo, dhi, hi)
    while 0 < f_hi < target and hi < 1e6:
        hi *= 2.0
        f_hi = _number_factor(dlo, dhi, hi)
    if f_hi <= 0 or target >= f_hi:
        return hi
    for _ in range(50):
        mid = math.sqrt(lo * hi)
        if _number_factor(dlo, dhi, mid) >= target:
            hi = mid
        else:
            lo = mid
    return math.sqrt(lo * hi)


def number_activation_fraction(dcrit, dlo, dhi, ratio):
    if dlo <= 0 or dhi <= dlo:
        return 0.0
    if dcrit <= dlo:
        return 1.0
    if dcrit >= dhi:
        return 0.0
    tot = _inv_d3_integral(dlo, dhi, ratio, dlo)
    act = _inv_d3_integral(dlo, dhi, ratio, dcrit)
    return max(0.0, min(1.0, act / tot)) if tot > 0 else 0.0


# ---------------------------------------------------------------- thermodynamics
def svp_water_goffgratch(T):
    """Goff & Gratch (1946) saturation vapour pressure over water (Pa), CAM default."""
    ts = 373.16
    e1 = 11.344 * (1.0 - T / ts)
    e2 = -3.49149 * (ts / T - 1.0)
    f1 = -7.90298 * (ts / T - 1.0)
    f2 = 5.02808 * math.log10(ts / T)
    f3 = -1.3816e-7 * (10.0 ** e1 - 1.0)
    f4 = 8.1328e-3 * (10.0 ** e2 - 1.0)
    f5 = math.log10(1013.246)
    return 10.0 ** (f1 + f2 + f3 + f4 + f5) * 100.0


def qsat_water(T, P):
    es = min(svp_water_goffgratch(T), P)
    return es, EPSILO * es / max(P - (1.0 - EPSILO) * es, 1e-30)


# ---------------------------------------------------------------- modes
def modes(pop: Populations, custom_kappa: dict | None = None) -> list[dict]:
    """Active (bin, class) modes with their critical supersaturation."""
    T, P = pop.T, pop.P
    rho_air = P / (RAIR * T)
    aten = 2.0 * MWH2O * SURFTEN / (R_UNIV * T * RHOH2O)
    kap = np.array([kappa_of(s, custom_kappa) for s in pop.species])
    iw = pop.iw
    rlo, rup = pop.r_low_up_cm()
    rmid = pop.r_mid_cm()
    out = []
    for b in range(pop.nb):
        r_m = rmid[b] * 1e-2
        for c in range(pop.nc):
            n = pop.number[b, c]
            if n <= 0 or n * 1e-6 <= NUMBER_FLOOR_CM3:
                continue
            m = np.maximum(pop.mass[b, c], 0.0)
            m[iw] = 0.0
            vol = m / pop.density_gcm3                          # ug/m3 / (g/cm3) = 1e-12 m3/m3
            vtot = vol.sum()
            if vtot <= 0:
                continue
            keff = float((vol * kap).sum() / vtot)
            if keff <= 0:
                continue
            sc = 2.0 * aten * math.sqrt(aten / (27.0 * keff * r_m ** 3))
            dryvol = vtot * 1e-12                               # m3/m3
            dlo, dhi = 2 * rlo[b] * 1e-2, 2 * rup[b] * 1e-2
            out.append({"bin": b, "cls": c, "n_m3": n, "kappa": keff, "sc": sc,
                        "dlo": dlo, "dhi": dhi, "dmid": 2 * r_m,
                        "tilt": tilt_ratio(dlo, dhi, n, dryvol)})
    return out


def _activated(md, s):
    if s <= 0:
        return 0.0
    dcrit = md["dmid"] * (md["sc"] / s) ** (2.0 / 3.0)
    return md["n_m3"] * number_activation_fraction(dcrit, md["dlo"], md["dhi"], md["tilt"])


def smax_arg(pop: Populations, md: list[dict], w: float) -> float | None:
    """Abdul-Razzak & Ghan maximum supersaturation (fraction) for updraft w (m/s)."""
    T, P = pop.T, pop.P
    rho_air = P / (RAIR * T)
    aten = 2.0 * MWH2O * SURFTEN / (R_UNIV * T * RHOH2O)
    w = max(w, WSUB_MIN)
    diff0 = 0.211e-4 * (PREF / P) * (T / TMELT) ** 1.94
    conduct0 = (5.69 + 0.017 * (T - TMELT)) * 4.186e2 * 1e-5
    _, qs = qsat_water(T, P)
    dqsdt = LATVAP / (RH2O * T * T) * qs
    alpha = GRAVIT * (LATVAP / (CPAIR * RH2O * T * T) - 1.0 / (RAIR * T))
    gamma = (1.0 + LATVAP / CPAIR * dqsdt) / (rho_air * qs)
    G = 1.0 / (RHOH2O / (diff0 * rho_air * qs) + LATVAP * RHOH2O / (conduct0 * T) * (LATVAP / (RH2O * T) - 1.0))
    beta = 4.0 * math.pi * RHOH2O * G * gamma
    if not (alpha > 0 and G > 0 and beta > 0):
        return None
    sqrtG = math.sqrt(G)
    alw = alpha * w
    sqrtalw = math.sqrt(alw)
    zeta = 2.0 * sqrtalw * aten / (3.0 * sqrtG)
    tot = 0.0
    for m in md:
        eta = 2.0 * alw * sqrtalw / (m["n_m3"] * beta * sqrtG) if m["n_m3"] > 1e-20 else 1e20
        g1 = math.sqrt(zeta / eta) ** 3
        g2 = math.sqrt(m["sc"] / math.sqrt(eta + 3.0 * zeta)) ** 3
        tot += (F1_BIN * g1 + F2_BIN * g2) / (m["sc"] * m["sc"])
    return 1.0 / math.sqrt(tot) if tot > 0 else None


def compute(pop: Populations, opts: dict | None = None, custom_kappa: dict | None = None) -> dict:
    o = dict(DEFAULTS, **(opts or {}))
    md = modes(pop, custom_kappa)
    ntot = float(pop.number.sum() * 1e-6)
    spec = []
    for s_pct in o["supersaturations_pct"]:
        s = s_pct / 100.0
        by_cls = np.zeros(pop.nc)
        for m in md:
            by_cls[m["cls"]] += _activated(m, s)
        spec.append({"S_pct": s_pct, "ccn_cm3": float(by_cls.sum() * 1e-6),
                     "by_class_cm3": (by_cls * 1e-6).tolist()})
    act = []
    for w in o["updrafts_m_s"]:
        smax = smax_arg(pop, md, w) if md else None
        by_cls = np.zeros(pop.nc)
        by_bin = np.zeros(pop.nb)
        if smax:
            for m in md:
                a = _activated(m, smax)
                by_cls[m["cls"]] += a
                by_bin[m["bin"]] += a
        n_cls = pop.number.sum(axis=0)
        act.append({"w_m_s": w, "smax_pct": (smax * 100 if smax else None),
                    "nd_cm3": float(by_cls.sum() * 1e-6),
                    "by_class_cm3": (by_cls * 1e-6).tolist(),
                    "actfrac_by_class": [float(a / n) if n > 0 else 0.0 for a, n in zip(by_cls, n_cls)],
                    "by_bin_cm3": (by_bin * 1e-6).tolist()})
    kappa_cls = []
    for c in range(pop.nc):
        ks = [(m["n_m3"], m["kappa"]) for m in md if m["cls"] == c]
        n = sum(k[0] for k in ks)
        kappa_cls.append(float(sum(k[0] * k[1] for k in ks) / n) if n > 0 else None)
    return {"n_total_cm3": ntot, "spectrum": spec, "activation": act,
            "kappa_by_class": kappa_cls, "options": o}


def series(run: dict, opts: dict | None = None, every: int = 1) -> dict:
    from .common import populations
    times = run["results"]["time_s"]
    idx = list(range(0, len(times), max(1, every)))
    if idx[-1] != len(times) - 1:
        idx.append(len(times) - 1)
    custom = {c["name"]: c.get("kappa") for c in run["case"]["species"].get("custom", []) if "kappa" in c}
    out = {"time_s": [float(times[i]) for i in idx], "resolved": [], "internal_mean": []}
    for i in idx:
        p = populations(run, i)
        out["resolved"].append(compute(p, opts, custom))
        out["internal_mean"].append(compute(p.internal_mean(), opts, custom))
    return out
