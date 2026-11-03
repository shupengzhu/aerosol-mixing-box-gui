#!/usr/bin/env python3
"""Adiabatic parcel reference for droplet activation (constant updraft).

Same thermodynamic coefficients as CAM ndrop/activate_modal (alpha, gamma, G), so the
parcel and the Abdul-Razzak & Ghan (ARG) parameterisation differ only in how the
supersaturation balance is solved:
  dS/dt   = alpha w - gamma * 4 pi rho_w * sum_i N_i r_i^2 dr_i/dt
  dr_i/dt = G / r_i * (S - S_eq,i(r_i)),  S_eq = exp(A/r) (r^3 - rd^3) / (r^3 - rd^3 (1 - kappa)) - 1
Particles start at Koehler equilibrium for S0 (default -0.5 %); the run stops at the
supersaturation maximum. Activated = number with critical supersaturation <= S_max
(the ARG definition).

parcel_nd(T, P, w, n_m3, rd_m, kappa) -> (smax, nd_m3); group_particles() reduces a
particle population to (log rd, kappa) groups first.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.sparse import lil_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scrambox.post.ccn import (CPAIR, GRAVIT, LATVAP, MWH2O, PREF, R_UNIV, RAIR,  # noqa: E402
                               RH2O, RHOH2O, SURFTEN, TMELT, qsat_water)

KAPPA_EDGES = np.array([0, 1e-6, 0.01, 0.03, 0.06, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0, 1.3])


def coefficients(T, P):
    rho_air = P / (RAIR * T)
    diff0 = 0.211e-4 * (PREF / P) * (T / TMELT) ** 1.94
    conduct0 = (5.69 + 0.017 * (T - TMELT)) * 4.186e2 * 1e-5
    _, qs = qsat_water(T, P)
    dqsdt = LATVAP / (RH2O * T * T) * qs
    alpha = GRAVIT * (LATVAP / (CPAIR * RH2O * T * T) - 1.0 / (RAIR * T))
    gamma = (1.0 + LATVAP / CPAIR * dqsdt) / (rho_air * qs)
    G = 1.0 / (RHOH2O / (diff0 * rho_air * qs) + LATVAP * RHOH2O / (conduct0 * T) * (LATVAP / (RH2O * T) - 1.0))
    A = 2.0 * MWH2O * SURFTEN / (R_UNIV * T * RHOH2O)
    return alpha, gamma, G, A


def s_crit(A, rd, kappa):
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(kappa > 1e-8, 2 * A * np.sqrt(A / (27 * kappa * rd ** 3)), np.inf)


def s_eq(r, rd, kappa, A):
    r3, d3 = r ** 3, rd ** 3
    return np.exp(A / r) * (r3 - d3) / (r3 - d3 * (1 - kappa)) - 1.0


def equilibrium_radius(S0, rd, kappa, A):
    out = rd.copy() * 1.0001
    for i in range(len(rd)):
        f = lambda r: s_eq(r, rd[i], kappa[i], A) - S0
        lo, hi = rd[i] * (1 + 1e-9), rd[i] * 1e3
        if kappa[i] > 1e-8 and f(lo) < 0 < f(hi):
            out[i] = brentq(f, lo, hi, xtol=rd[i] * 1e-9)
    return out


def group_particles(rd_m, kappa, n_m3, nsize=120):
    """(log rd, kappa) groups: number, volume-mean dry radius, volume-weighted kappa."""
    edges = np.logspace(-9.5, -4.5, nsize + 1)
    ib = np.clip(np.searchsorted(edges, rd_m) - 1, 0, nsize - 1)
    ik = np.clip(np.searchsorted(KAPPA_EDGES, kappa, side="right") - 1, 0, len(KAPPA_EDGES) - 2)
    key = ib * 100 + ik
    out_n, out_r, out_k = [], [], []
    for k in np.unique(key):
        s = key == k
        N = n_m3[s].sum()
        v = (n_m3[s] * rd_m[s] ** 3).sum()
        out_n.append(N)
        out_r.append((v / N) ** (1 / 3))
        out_k.append((n_m3[s] * rd_m[s] ** 3 * kappa[s]).sum() / v)
    return np.array(out_n), np.array(out_r), np.array(out_k)


def parcel_nd(T, P, w, n_m3, rd_m, kappa, S0=-0.005):
    alpha, gamma, G, A = coefficients(T, P)
    keep = n_m3 > 0
    n, rd, kap = n_m3[keep], rd_m[keep], np.maximum(kappa[keep], 0.0)
    r0 = equilibrium_radius(S0, rd, kap, A)
    m = len(n)
    coef = gamma * 4 * math.pi * RHOH2O

    def rhs(t, y):
        S, r = y[0], np.maximum(y[1:], rd * (1 + 1e-12))
        drdt = G / r * (S - s_eq(r, rd, kap, A))
        dS = alpha * w - coef * np.sum(n * r * r * drdt)
        return np.concatenate(([dS], drdt))

    sp = lil_matrix((m + 1, m + 1))
    sp[0, :] = 1
    sp[:, 0] = 1
    sp.setdiag(1)

    def peak(t, y):
        return rhs(t, y)[0]
    peak.terminal, peak.direction = True, -1

    t_end = 20 * (abs(S0) + 0.02) / (alpha * w)
    sol = solve_ivp(rhs, (0, t_end), np.concatenate(([S0], r0)), method="BDF", jac_sparsity=sp.tocsr(),
                    events=peak, rtol=1e-6, atol=np.concatenate(([1e-9], rd * 1e-6)), max_step=t_end / 50)
    S = sol.y[0]
    smax = float(S.max())
    nd = float(n[s_crit(A, rd, kap) <= smax].sum())
    return smax, nd


if __name__ == "__main__":
    # self-test: one lognormal mode, parcel vs modal ARG and the same mode as k monodisperse sections
    from scipy.special import erfc
    T, P, w = 283.15, 9.0e4, 0.5
    N, dg, sg, kap = 1000e6, 0.1e-6, 1.6, 0.5
    alpha, gamma, G, A = coefficients(T, P)
    lsg = math.log(sg)
    edges = np.linspace(math.log(dg) - 5 * lsg, math.log(dg) + 5 * lsg, 401)
    mid = 0.5 * (edges[1:] + edges[:-1])
    frac = 0.5 * (erfc((edges[:-1] - math.log(dg)) / (math.sqrt(2) * lsg)) - erfc((edges[1:] - math.log(dg)) / (math.sqrt(2) * lsg)))
    smax, nd = parcel_nd(T, P, w, N * frac, np.exp(mid) / 2, np.full(400, kap))
    print(f"parcel: smax {smax * 100:.4f} %  Nd {nd / 1e6:.1f} cm-3")
