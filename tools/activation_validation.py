#!/usr/bin/env python3
"""Validation of S_max schemes for SCRAM on the production 7-bin grid against the parcel model.

Synthetic aerosol: unimodal, bimodal (Aitken + accumulation, different kappa) and external
mixtures (hygroscopic + nearly hydrophobic particles of the same size). Each component is one
composition class. Its lognormal is put on the production bins (number and dry volume per
bin x class, as SCRAM stores it); activation after S_max uses SCRAM's rule (critical
supersaturation at the bin mid radius, linear dry-mass tilt inside the bin).

  truth     parcel model on the continuous lognormals (400 sections per component)
  modal     ARG 2000 on the true lognormal modes (what MAM4 would do with perfect modes)
  current   modal maxsat over every bin x class slot, f1 0.5, f2 1 (scram_cloud_aerosol_export)
  B         one effective population (ARG 2002 effective critical supersaturation, ln sigma
            from the spread of ln sc)
  C<w>      slots grouped by a fixed grid in log10(sc) of width w decades; each group is one
            effective mode (as B); modal maxsat over the groups
  B/C/current are also evaluated with every class duplicated into two identical half classes
  (split test: must not change).

Usage: python3 tools/activation_validation.py [out_dir]   (uses all cores)
"""
from __future__ import annotations

import itertools
import json
import math
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parcel_activation import coefficients, parcel_nd                     # noqa: E402
from scrambox.post.ccn import number_activation_fraction, tilt_ratio      # noqa: E402

BOUNDS_UM = np.array([0.003, 0.01, 0.0398, 0.1585, 0.6310, 2.5, 5.0, 10.0])   # production cfg
T, P = 283.15, 9.0e4
SQ2 = math.sqrt(2)
UPDRAFTS = (0.1, 0.3, 1.0, 3.0)
C_WIDTHS = (1.0,)
SCHEMES = ["current", "Bt", "FN", "FNr"]


def cdf(x_ln_d, dg, sg):
    return 0.5 * erfc(-(x_ln_d - math.log(dg)) / (SQ2 * math.log(sg)))


def slots(components):
    """components: list of (N m-3, dg m, sg, kappa) -> per bin x class N, V, dlo, dhi, kappa, cls."""
    out = []
    e = np.log(BOUNDS_UM * 1e-6)
    for c, (N, dg, sg, kap) in enumerate(components):
        lsg = math.log(sg)
        for b in range(len(BOUNDS_UM) - 1):
            n = N * (cdf(e[b + 1], dg, sg) - cdf(e[b], dg, sg))
            if n < 1e3:
                continue
            # volume in the bin: third-moment lognormal (dg_v = dg exp(3 ln^2 sg))
            dgv = dg * math.exp(3 * lsg ** 2)
            v = N * math.pi / 6 * dg ** 3 * math.exp(4.5 * lsg ** 2) * (cdf(e[b + 1], dgv, sg) - cdf(e[b], dgv, sg))
            out.append({"N": n, "V": v, "dlo": BOUNDS_UM[b] * 1e-6, "dhi": BOUNDS_UM[b + 1] * 1e-6, "kappa": kap, "cls": c})
    return out


def prepare(sl, A):
    for s in sl:
        s["dmid"] = math.sqrt(s["dlo"] * s["dhi"])
        s["sc"] = 2 * A * math.sqrt(A / (27 * s["kappa"] * (s["dmid"] / 2) ** 3))
        s["tilt"] = tilt_ratio(s["dlo"], s["dhi"], s["N"], s["V"])
    return sl


def activated(sl, smax):
    tot = 0.0
    for s in sl:
        dcrit = s["dmid"] * (s["sc"] / smax) ** (2 / 3)
        tot += s["N"] * number_activation_fraction(dcrit, s["dlo"], s["dhi"], s["tilt"])
    return tot


def arg_consts(w):
    alpha, gamma, G, A = coefficients(T, P)
    beta = 2 * math.pi * 1000.0 * G * gamma
    alw = alpha * w
    return 2 / 3 * math.sqrt(alw) * A / math.sqrt(G), alw ** 1.5 / (beta * math.sqrt(G)), A


def maxsat(N, sm, f1, f2, zeta, pref):
    N, sm, f1, f2 = map(np.asarray, (N, sm, f1, f2))
    eta = pref / N
    tot = np.sum((f1 * (zeta / eta) ** 1.5 + f2 * (sm ** 2 / (eta + 3 * zeta)) ** 0.75) / sm ** 2)
    return 1 / math.sqrt(tot)


def eff_mode(group, geo=False):
    """Effective lognormal mode of a group of slots. geo=False: ARG-2002 effective critical
    supersaturation (N / sum N sc^-2/3)^1.5 (number-mean radius); geo=True: number-weighted
    geometric mean of sc (median radius, as the modal formula expects)."""
    N = np.array([s["N"] for s in group])
    sc = np.array([s["sc"] for s in group])
    Nt = N.sum()
    ls = np.log(sc)
    sm = math.exp(np.sum(N * ls) / Nt) if geo else (Nt / np.sum(N * sc ** (-2 / 3))) ** 1.5
    var = np.sum(N * ls ** 2) / Nt - (np.sum(N * ls) / Nt) ** 2
    # within-bin spread of ln sc for a log-uniform bin: (1.5 ln(dhi/dlo))^2 / 12
    var += np.sum(N * np.array([(1.5 * math.log(s["dhi"] / s["dlo"])) ** 2 / 12 for s in group])) / Nt
    lnsig = math.sqrt(max(var, 1e-4)) / 1.5
    return Nt, sm, lnsig


def tilt_moments(s, nsub=24):
    """Mean and variance of ln d inside a bin under SCRAM's linear dry-mass tilt."""
    e = np.geomspace(s["dlo"], s["dhi"], nsub + 1)
    F = np.array([number_activation_fraction(x, s["dlo"], s["dhi"], s["tilt"]) for x in e])
    wgt = np.maximum(F[:-1] - F[1:], 0.0)
    x = np.log(np.sqrt(e[:-1] * e[1:]))
    wgt = wgt / wgt.sum()
    m = float(np.sum(wgt * x))
    return m, float(np.sum(wgt * (x - m) ** 2)) + (math.log(e[1] / e[0])) ** 2 / 12


def eff_mode_tilt(group, A):
    """One effective lognormal from the exact first and second moments of ln sc of all slots
    (ln sc = ln(2A/3 sqrt(4A/(3 kappa))) - 1.5 ln d); every sum is additive over slots."""
    Nt = sum(s["N"] for s in group)
    m1 = m2 = 0.0
    for s in group:
        mu, var = tilt_moments(s)
        c = math.log(2 * A * math.sqrt(A / (27 * s["kappa"]))) + 1.5 * math.log(2.0)   # sc = c * d^-1.5
        mls = c - 1.5 * mu
        m1 += s["N"] * mls
        m2 += s["N"] * (mls ** 2 + 2.25 * var)
    mean = m1 / Nt
    lnsig = math.sqrt(max(m2 / Nt - mean ** 2, 1e-4)) / 1.5
    return Nt, math.exp(mean), lnsig


def fn_spectrum(sl, nsub=16):
    """Critical-supersaturation spectrum of all slots as discrete sub-sections (SCRAM tilt)."""
    N, sc = [], []
    for s in sl:
        e = np.geomspace(s["dlo"], s["dhi"], nsub + 1)
        F = np.array([number_activation_fraction(x, s["dlo"], s["dhi"], s["tilt"]) for x in e])
        N.append(s["N"] * np.maximum(F[:-1] - F[1:], 0.0))
        sc.append(s["sc"] * (s["dmid"] / np.sqrt(e[:-1] * e[1:])) ** 1.5)
    return np.concatenate(N), np.concatenate(sc)


def smax_fn(N, sc, w, T_, P_, a_emp=2.0, stats=None):
    """Nenes & Seinfeld (2003) / Fountoukis & Nenes (2005) water-vapour balance on a sectional
    CCN spectrum, with the parcel/CAM coefficients (r dr/dt = G S, Kelvin A in radius form):
      alpha w = gamma 4 pi rho_w G smax sum_i N_i r_i,
      r_i = sqrt(G (smax^2 - sc_i^2) / (alpha w))  for sc_i <= s_part (grown CCN)
      r_i = 2 A / (3 sc_i)                          for s_part < sc_i <= smax (critical radius)
    s_part from the discriminant D = smax^4 - 16 A^2 alpha w / (9 G) (FN05 sec. 3.4); D < 0:
    s_part = smax min(2e7 (2A) / 3 smax^-0.3824, 1) with the diameter-form A of FN05."""
    alpha, gamma, G, A = coefficients(T_, P_)
    aw = alpha * w
    coef = gamma * 4 * math.pi * 1000.0 * G

    def f(smax):
        D = smax ** 4 - 16 * A ** 2 * aw / (9 * G)
        if D >= 0:
            spart = smax * math.sqrt(0.5 * (1 + math.sqrt(1 - 16 * A ** 2 * aw / (9 * G * smax ** 4))))
        else:
            spart = smax * min(2e7 * (a_emp * A) / 3 * smax ** -0.3824, 1.0)
        m1 = sc <= spart
        m2 = (sc > spart) & (sc <= smax)
        I = np.sum(N[m1] * np.sqrt(G * (smax ** 2 - sc[m1] ** 2) / aw)) + np.sum(N[m2] * 2 * A / (3 * sc[m2]))
        return aw - coef * smax * I

    lo, hi = 1e-7, 0.5
    if f(hi) > 0:
        return hi
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    s = math.sqrt(lo * hi)
    if stats is not None:
        stats.append(s ** 4 - 16 * A ** 2 * aw / (9 * G) >= 0)
    return s


def smax_mbn(N, sc, w, T_, P_, a_emp=2.0, stats=None):
    """Morales Betancourt & Nenes (2014) revised population splitting on the sectional CCN
    spectrum: radius at s_max = D_pc/sqrt(3) (inertially limited, sc < s-), growth term
    (s- <= sc <= s+) or critical radius (s+ < sc <= s_max); for s_max < xi_c both roots from
    their Eq. (14), continuous with Eq. (10) at s_max = xi_c. Same coefficients as smax_fn."""
    alpha, gamma, G, A = coefficients(T_, P_)
    aw = alpha * w
    coef = gamma * 4 * math.pi * 1000.0 * G
    xi = (16 * A ** 2 * aw / (9 * G)) ** 0.25
    Ad = a_emp * A
    rc = 2 * A / (3 * sc)

    def split(smax):
        if smax > xi:
            q = math.sqrt(1 - (xi / smax) ** 4)
            return smax * math.sqrt(0.5 * (1 + q)), smax * math.sqrt(0.5 * (1 - q))
        r = min(2e7 * Ad / 3 * (smax ** -0.3824 - xi ** -0.3824) + 1 / math.sqrt(2), 1.0)
        return smax * r, smax * r

    def f(smax):
        sp, sm = split(smax)
        act = sc <= smax
        rg = np.sqrt(np.maximum(G * (smax ** 2 - sc ** 2) / aw, 0.0))
        ri = np.where(sc < sm, rc / math.sqrt(3), np.where(sc <= sp, rg, rc))
        return aw - coef * smax * np.sum(N[act] * ri[act])

    lo, hi = 1e-7, 0.5
    if f(hi) > 0:
        return hi
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if f(mid) > 0:
            lo = mid
        else:
            hi = mid
    s = math.sqrt(lo * hi)
    if stats is not None:
        stats.append(s > xi)
    return s


def smax_scheme(name, sl, zeta, pref, A=None):
    if name == "MBN":
        N, sc = fn_spectrum(sl)
        return smax_mbn(N, sc, smax_scheme.w, T, P)
    if name in ("FN", "FNr"):
        N, sc = fn_spectrum(sl)
        return smax_fn(N, sc, smax_scheme.w, T, P, a_emp=2.0 if name == "FN" else 1.0)
    if name == "current":
        return maxsat([s["N"] for s in sl], [s["sc"] for s in sl], 0.5, 1.0, zeta, pref)
    if name == "Bt":
        Nt, sm, lns = eff_mode_tilt(sl, A)
        return maxsat([Nt], [sm], 0.5 * math.exp(2.5 * lns ** 2), 1 + 0.25 * lns, zeta, pref)
    geo = name.endswith("g")
    base = name[:-1] if geo else name
    if base == "B":
        groups = [sl]
    else:
        width = float(base[1:])
        key = lambda s: math.floor(math.log10(s["sc"]) / width)
        groups = [list(g) for _, g in itertools.groupby(sorted(sl, key=key), key=key)]
    modes = [eff_mode(g, geo) for g in groups]
    lns = np.array([m[2] for m in modes])
    return maxsat([m[0] for m in modes], [m[1] for m in modes], 0.5 * np.exp(2.5 * lns ** 2), 1 + 0.25 * lns, zeta, pref)


def case_eval(args):
    label, comps, w = args
    zeta, pref, A = arg_consts(w)
    smax_scheme.w = w
    # truth: parcel on fine sections of every component
    ns, rs, ks = [], [], []
    for N, dg, sg, kap in comps:
        lsg = math.log(sg)
        e = np.linspace(math.log(dg) - 5 * lsg, math.log(dg) + 5 * lsg, 401)
        ns.append(N * np.diff(cdf(e, dg, sg))); rs.append(np.exp(0.5 * (e[1:] + e[:-1])) / 2); ks.append(np.full(400, kap))
    s_true, nd_true = parcel_nd(T, P, w, np.concatenate(ns), np.concatenate(rs), np.concatenate(ks))
    # modal ARG on the true modes
    sm = [2 * A * math.sqrt(A / (27 * k * (dg / 2) ** 3)) for N, dg, sg, k in comps]
    lns = np.array([math.log(sg) for _, _, sg, _ in comps])
    s_mod = maxsat([c[0] for c in comps], sm, 0.5 * np.exp(2.5 * lns ** 2), 1 + 0.25 * lns, zeta, pref)
    nd_mod = sum(N * 0.5 * erfc(2 / 3 * math.log(smi / s_mod) / (SQ2 * math.log(sg))) for (N, dg, sg, k), smi in zip(comps, sm))
    sl = prepare(slots(comps), A)
    split = []
    for s in sl:
        for h in range(2):
            split.append(dict(s, N=s["N"] / 2, V=s["V"] / 2, cls=2 * s["cls"] + h))
    # SCRAM's counting rule with the exact (parcel) S_max: error of the 7-bin representation alone
    res = {"label": label, "w": w, "truth": nd_true / 1e6, "modal": nd_mod / 1e6, "smax_true": s_true,
           "smax_modal": s_mod, "count_at_true_smax": activated(sl, s_true) / 1e6}
    for name in SCHEMES:
        s = smax_scheme(name, sl, zeta, pref, A)
        res[name] = activated(sl, s) / 1e6
        res["smax_" + name] = s
        res[name + "_split"] = activated(split, smax_scheme(name, split, zeta, pref, A)) / 1e6
    return res


def cases():
    out = []
    for N, dg, sg, k in itertools.product((50e6, 200e6, 1000e6, 5000e6), (0.04e-6, 0.08e-6, 0.15e-6), (1.5, 1.8, 2.2), (0.1, 0.3, 0.6, 1.2)):
        out.append((f"uni N{N / 1e6:.0f} dg{dg * 1e6:.2f} sg{sg} k{k}", [(N, dg, sg, k)]))
    for Na, Nc, (ka, kc) in itertools.product((200e6, 1000e6, 3000e6), (50e6, 300e6, 1000e6), ((0.2, 0.5), (0.6, 0.6), (0.1, 1.0))):
        out.append((f"bi Nait{Na / 1e6:.0f} Nacc{Nc / 1e6:.0f} k{ka}/{kc}", [(Na, 0.03e-6, 1.6, ka), (Nc, 0.15e-6, 1.6, kc)]))
    for N, f, dg in itertools.product((200e6, 1000e6, 3000e6), (0.3, 0.7), (0.08e-6, 0.15e-6)):
        out.append((f"ext N{N / 1e6:.0f} hydrophobic{f} dg{dg * 1e6:.2f}", [(N * (1 - f), dg, 1.6, 0.5), (N * f, dg, 1.6, 0.01)]))
    return [(lab, comps, w) for lab, comps in out for w in UPDRAFTS]


def summarize(res):
    names = ["count_at_true_smax", "modal"] + SCHEMES
    rr = [r for r in res if r["truth"] > 1]
    L0 = ["## S_max ratio scheme / parcel (all cases)", "", "| scheme | median | 10-90 % |", "|---|---|---|"]
    for n in ["modal"] + SCHEMES:
        x = np.array([r["smax_" + n] / r["smax_true"] for r in rr])
        L0.append(f"| {n} | {np.median(x):.2f} | {np.percentile(x, 10):.2f}-{np.percentile(x, 90):.2f} |")
    L0.append("")
    L = ["# SCRAM S_max schemes on the production 7-bin grid vs parcel", "",
         f"{len(res)} cases (unimodal, bimodal, external mixtures; w {UPDRAFTS} m/s). Ratio = scheme Nd / parcel Nd "
         "(cases with parcel Nd > 1 cm-3).", ""]
    for sub, sel in (("all", lambda r: True), ("unimodal", lambda r: r["label"].startswith("uni")),
                     ("bimodal", lambda r: r["label"].startswith("bi")), ("external", lambda r: r["label"].startswith("ext")),
                     ("N <= 1000 cm-3, w >= 0.3", lambda r: ("N5000" not in r["label"] and "3000" not in r["label"]) and r["w"] >= 0.3)):
        rr = [r for r in res if sel(r) and r["truth"] > 1]
        L += [f"## {sub} ({len(rr)} cases)", "", "| scheme | median ratio | 10-90 % | RMS log10 error | split change (median, max) |",
              "|---|---|---|---|---|"]
        for n in names:
            x = np.array([r[n] / r["truth"] for r in rr])
            x = np.clip(x, 1e-3, None)
            sp = ("-" if n in ("modal", "count_at_true_smax") else
                  f"{np.median([abs(r[n + '_split'] / max(r[n], 1e-9) - 1) for r in rr]) * 100:.0f}%, "
                  f"{np.max([abs(r[n + '_split'] / max(r[n], 1e-9) - 1) for r in rr]) * 100:.0f}%")
            L.append(f"| {n} | {np.median(x):.2f} | {np.percentile(x, 10):.2f}-{np.percentile(x, 90):.2f} | "
                     f"{math.sqrt(np.mean(np.log10(x) ** 2)):.3f} | {sp} |")
        L.append("")
    return L + L0


def main():
    out = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    cs = cases()
    with Pool(min(96, len(cs))) as p:
        res = p.map(case_eval, cs, chunksize=2)
    (out / "activation_validation.json").write_text(json.dumps(res))
    L = summarize(res)
    (out / "activation_validation.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
