#!/usr/bin/env python3
"""Candidate S_max schemes for SCRAM's bin x class populations, tested against the parcel model.

  current  modal ARG maxsat over every bin x class slot, f1 0.5, f2 1 (scram_cloud_aerosol_export)
  A_types  WRF-Chem module_mixactivate iquasisect_option 3/4 with types = composition classes:
           per class, all bins form one mode (N_c, ln sigma from the number-weighted spread of
           ln r, ARG-2002 effective critical supersaturation
           sm_c = (N_c / sum_i N_i sm_i^(-2/3))^(3/2)); modal maxsat over the classes
  B_single the same with all bins of all classes as one population; ln sigma from the
           number-weighted spread of ln sm (= 1.5 ln sigma_r for uniform kappa), so every sum
           is additive and the result cannot depend on how populations are split
Activation after S_max is unchanged: every bin x class slot activates with its own critical
supersaturation and SCRAM's sub-bin tilt (scrambox.post.ccn._activated).

Tests: (1) one lognormal mode cut into k sections; (2) identical populations split into k
copies; (3) WP2 end-of-run states (SCRAM box runs and PartMC->NC1/NC6/NC11) against the
parcel truth from runs/wp2/activation_parcel.json.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from parcel_activation import coefficients, parcel_nd                 # noqa: E402
from scrambox.post import ccn as C                                     # noqa: E402
from scrambox.post.common import populations                           # noqa: E402
from scrambox.run import load_run                                      # noqa: E402

SQ2 = math.sqrt(2)
LNSIG_MIN = 0.01                     # WRF: gmlnsig = sqrt(max(1e-4, var))


def _arg(T, P, w):
    alpha, gamma, G, A = coefficients(T, P)
    beta = 2 * math.pi * 1000.0 * G * gamma
    alw = alpha * w
    zeta = 2 / 3 * math.sqrt(alw) * A / math.sqrt(G)
    return zeta, alw ** 1.5 / (beta * math.sqrt(G))


def maxsat(N, sm, lnsig, zeta, pref, f1=None, f2=None):
    N, sm, lnsig = map(np.atleast_1d, (N, sm, lnsig))
    ok = N > 1e-20
    eta = pref / N[ok]
    f1 = 0.5 * np.exp(2.5 * lnsig[ok] ** 2) if f1 is None else np.broadcast_to(f1, N.shape)[ok]
    f2 = 1 + 0.25 * lnsig[ok] if f2 is None else np.broadcast_to(f2, N.shape)[ok]
    tot = np.sum((f1 * (zeta / eta) ** 1.5 + f2 * (sm[ok] ** 2 / (eta + 3 * zeta)) ** 0.75) / sm[ok] ** 2)
    return 1 / math.sqrt(tot)


def group_stats(N, sm, r):
    """ARG-2002 effective critical supersaturation and two ln sigma estimates for a group."""
    Nt = N.sum()
    sm_eff = (Nt / np.sum(N * sm ** (-2 / 3))) ** 1.5
    lr = np.log(r)
    lnsig_r = math.sqrt(max(LNSIG_MIN ** 2, np.sum(N * lr ** 2) / Nt - (np.sum(N * lr) / Nt) ** 2))
    ls = np.log(sm)
    lnsig_s = math.sqrt(max((1.5 * LNSIG_MIN) ** 2, np.sum(N * ls ** 2) / Nt - (np.sum(N * ls) / Nt) ** 2)) / 1.5
    return Nt, sm_eff, lnsig_r, lnsig_s


def smax_schemes(N, sm, r, cls, T, P, w):
    """N (m-3), sm (critical supersaturation at the slot's mean radius), r (m), cls (class id)."""
    zeta, pref = _arg(T, P, w)
    out = {"current": maxsat(N, sm, np.zeros_like(N), zeta, pref, 0.5, 1.0)}
    g = [group_stats(N[cls == c], sm[cls == c], r[cls == c]) for c in np.unique(cls)]
    out["A_types"] = maxsat(np.array([x[0] for x in g]), np.array([x[1] for x in g]),
                            np.array([x[2] for x in g]), zeta, pref)
    Nt, sm_eff, _, lnsig_s = group_stats(N, sm, r)
    out["B_single"] = maxsat(np.array([Nt]), np.array([sm_eff]), np.array([lnsig_s]), zeta, pref)
    return out


# ---------------------------------------------------------------- test 1 and 2: splitting
def lognormal_sections(N, dg, sg, kap, k, T):
    _, _, _, A = coefficients(T, 1e5)
    lsg = math.log(sg)
    e = np.linspace(math.log(dg) - 3 * lsg, math.log(dg) + 3 * lsg, k + 1)
    fr = np.diff(0.5 * erfc(-(e - math.log(dg)) / (SQ2 * lsg)))
    d = np.exp(0.5 * (e[1:] + e[:-1]))
    sm = 2 * A * np.sqrt(A / (27 * kap * (d / 2) ** 3))
    return N * fr, sm, d / 2, A


def lognormal_nd(N, dg, sg, kap, smax, A):
    sm = 2 * A * math.sqrt(A / (27 * kap * (dg / 2) ** 3))
    return N * 0.5 * erfc(2 / 3 * math.log(sm / smax) / (SQ2 * math.log(sg)))


def split_tests():
    T, P = 283.15, 9.0e4
    L = ["## 1. One lognormal mode cut into k sections (N 1000 cm⁻³, dg 0.1 µm, σg 1.6, κ 0.5)", "",
         "| w | parcel | modal ARG | scheme | k=7 | k=14 | k=42 | k=200 |", "|---|---|---|---|---|---|---|---|"]
    N, dg, sg, kap = 1000e6, 0.1e-6, 1.6, 0.5
    for w in (0.2, 0.5, 1.0):
        Nf, smf, rf, A = lognormal_sections(N, dg, sg, kap, 400, T)
        _, ndp = parcel_nd(T, P, w, Nf, rf, np.full(len(Nf), kap))
        zeta, pref = _arg(T, P, w)
        sm0 = 2 * A * math.sqrt(A / (27 * kap * (dg / 2) ** 3))
        modal = lognormal_nd(N, dg, sg, kap, maxsat(N, sm0, math.log(sg), zeta, pref), A)
        rows = {}
        for k in (7, 14, 42, 200):
            n, sm, r, _ = lognormal_sections(N, dg, sg, kap, k, T)
            for name, s in smax_schemes(n, sm, r, np.zeros(k), T, P, w).items():
                rows.setdefault(name, []).append(lognormal_nd(N, dg, sg, kap, s, A))
        for name, v in rows.items():
            L.append(f"| {w} | {ndp / 1e6:.0f} | {modal / 1e6:.0f} | {name} | " + " | ".join(f"{x / 1e6:.0f}" for x in v) + " |")
    L += ["", "## 2. External mixture 500 κ0.5 + 500 κ0.01 (0.1 µm), each class split into k identical copies, w 0.5", "",
          "| scheme | k=1 | k=3 | k=21 |", "|---|---|---|---|"]
    w = 0.5
    _, _, _, A = coefficients(T, P)
    rows = {}
    for k in (1, 3, 21):
        n = np.full(2 * k, 500e6 / k)
        kap2 = np.concatenate([np.full(k, 0.5), np.full(k, 0.01)])
        r = np.full(2 * k, 0.05e-6)
        sm = 2 * A * np.sqrt(A / (27 * kap2 * r ** 3))
        cls = np.concatenate([np.zeros(k), np.ones(k)]) if True else None
        for scheme_cls, tag in ((cls, ""), (np.arange(2 * k), " (copies as classes)")):
            for name, s in smax_schemes(n, sm, r, scheme_cls, T, P, w).items():
                rows.setdefault(name + tag, []).append(f"{s * 100:.3f}% / {n[sm <= s].sum() / 1e6:.0f}")
    _, ndp = parcel_nd(T, P, w, np.array([500e6, 500e6]), np.full(2, 0.05e-6), np.array([0.5, 0.01]))
    for name, v in rows.items():
        L.append(f"| {name} | " + " | ".join(v) + " |")
    L.append(f"| parcel | {ndp / 1e6:.0f} (any k) | | |")
    return L


# ---------------------------------------------------------------- test 3: WP2 states
def state_nd(run_dir, w):
    pop = populations(load_run(run_dir), -1)
    md = C.modes(pop)
    N = np.array([m["n_m3"] for m in md])
    sm = np.array([m["sc"] for m in md])
    r = np.array([m["dmid"] / 2 for m in md])
    cls = np.array([m["cls"] for m in md])
    res = {}
    for name, s in smax_schemes(N, sm, r, cls, pop.T, pop.P, w).items():
        res[name] = sum(C._activated(m, s) for m in md) / 1e6
    return res


def wp2_tests(root):
    truth = json.loads((root / "activation_parcel.json").read_text())
    L = ["", "## 3. WP2 end-of-run states vs parcel truth (Nd cm⁻³, bias in brackets)", ""]
    for scn, res in truth.items():
        L += [f"### {scn}", "", "| w | truth | state | current | A_types | B_single |", "|---|---|---|---|---|---|"]
        for w, row in res.items():
            t = float(np.mean([x["nd"] for x in row["truth"]]))
            states = [("SCRAM NC1", root / scn / "scram_NC1"), ("SCRAM NC6", root / scn / "scram_NC6"),
                      ("SCRAM NC11", root / scn / "scram_NC11")]
            states += [(f"PartMC→{nc}", [root / scn / f"partmc_s{s}_{nc}" for s in (1, 2, 3)]) for nc in ("NC1", "NC6", "NC11")]
            for label, d in states:
                ds = d if isinstance(d, list) else [d]
                vals = [state_nd(x, float(w)) for x in ds if x.is_dir()]
                if not vals:
                    continue
                cells = []
                for name in ("current", "A_types", "B_single"):
                    v = float(np.mean([x[name] for x in vals]))
                    cells.append(f"{v:.4g} ({100 * (v / t - 1):+.0f}%)")
                L.append(f"| {w} | {t:.4g} | {label} | " + " | ".join(cells) + " |")
        L.append("")
    return L


def main():
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/wp2")
    L = ["# S_max scheme candidates for SCRAM activation", ""] + split_tests() + wp2_tests(root)
    (root / "arg_sectional_candidates.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
