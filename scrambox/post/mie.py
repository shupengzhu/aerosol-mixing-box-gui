"""Homogeneous-sphere Mie efficiencies (Bohren & Huffman 1983, BHMIE), vectorised over
particles. Used for the Optics_mixing 0 (per-type tables) and 1 (volume-mixed refractive
index) paths, which the host model evaluates with CARMA's Mie routine for spheres."""
from __future__ import annotations

import numpy as np


def mie_q(x, m):
    """Qext, Qsca, g for size parameters x (N,) and complex refractive indices m (N,)."""
    x = np.atleast_1d(np.asarray(x, dtype=float))
    m = np.atleast_1d(np.asarray(m, dtype=complex)) * np.ones_like(x)
    qext = np.zeros_like(x)
    qsca = np.zeros_like(x)
    gsca = np.zeros_like(x)
    ok = x > 1e-8
    if not ok.any():
        return qext, qsca, gsca
    x, m = x[ok], m[ok]
    y = m * x
    xstop = x + 4.0 * x ** (1.0 / 3.0) + 2.0
    nstop = int(np.ceil(xstop.max()))
    nmx = int(np.ceil(max(xstop.max(), np.abs(y).max()))) + 15
    # logarithmic derivative D_n(mx) by downward recurrence
    d = np.zeros((nmx + 1, x.size), dtype=complex)
    for n in range(nmx, 0, -1):
        en = n / y
        d[n - 1] = en - 1.0 / (d[n] + en)
    psi0, psi1 = np.cos(x), np.sin(x)
    chi0, chi1 = -np.sin(x), np.cos(x)
    xi1 = psi1 - 1j * chi1
    an1 = bn1 = np.zeros_like(y)
    qs = np.zeros_like(x)
    gs = np.zeros_like(x)
    qe = np.zeros_like(x)
    for n in range(1, nstop + 1):
        # each particle stops at its own nstop; beyond it the upward recurrences of
        # small particles overflow, so their terms are masked out (and errors silenced)
        live = n <= xstop
        fn = (2.0 * n + 1.0) / (n * (n + 1.0))
        with np.errstate(all="ignore"):
            psi = np.where(live, (2.0 * n - 1.0) * psi1 / x - psi0, 0.0)
            chi = np.where(live, (2.0 * n - 1.0) * chi1 / x - chi0, 1.0)
            xi = psi - 1j * chi
            dn = d[n]
            an = np.where(live, ((dn / m + n / x) * psi - psi1) / ((dn / m + n / x) * xi - xi1), 0.0)
            bn = np.where(live, ((m * dn + n / x) * psi - psi1) / ((m * dn + n / x) * xi - xi1), 0.0)
        qs += (2.0 * n + 1.0) * (np.abs(an) ** 2 + np.abs(bn) ** 2)
        qe += (2.0 * n + 1.0) * (an.real + bn.real)
        gs += fn * (an.real * bn.real + an.imag * bn.imag)
        if n > 1:
            gs += ((n - 1.0) * (n + 1.0) / n) * (an1.real * an.real + an1.imag * an.imag +
                                                 bn1.real * bn.real + bn1.imag * bn.imag)
        psi0, psi1 = psi1, psi
        chi0, chi1 = chi1, chi
        xi1 = psi1 - 1j * chi1
        an1, bn1 = an, bn
    qsca_ok = 2.0 * qs / x ** 2
    qext_ok = 2.0 * qe / x ** 2
    g_ok = np.where(qs > 0, 2.0 * gs / np.where(qs > 0, qs, 1.0), 0.0)
    qext[ok], qsca[ok], gsca[ok] = qext_ok, np.minimum(qsca_ok, qext_ok), g_ok
    return qext, qsca, gsca
