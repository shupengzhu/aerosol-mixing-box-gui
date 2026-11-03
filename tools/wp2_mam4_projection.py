#!/usr/bin/env python3
"""WP2: how far is the MAM4 modal representation from PartMC, compared with SCRAM?

PartMC particles are mapped onto the four CESM2 MAM4 modes with CAM's own rules and the
CCN spectrum and droplet number are computed with CAM's formulas (ndrop ccncalc and
activate_modal/maxsat, single updraft). The same particles also give the particle-resolved
truth (per-particle kappa-Koehler CCN; ARG with every particle as its own population). This
is a static representation test: MAM4 receives the exact PartMC aerosol, so the result is
the error of the modal representation alone, not of MAM4's process code.

Mode assignment per particle (dry diameter d, volumes v_i, MAM4 physprops kappa):
  coarse   d >= 1 um (MAM4 emits dust and sea salt above 1 um into the coarse mode)
  pcarbon  BC or POA present, no dust/salt, d < 1 um and the SO4-equivalent coating
           sum(v_i kappa_i / kappa_SO4) below pi d_core^2 * dr, dr = 8 monolayers * 0.476 nm
           (n_so4_monolayers_pcage, modal_aero_coag/gasaerexch ageing rule)
  aitken   d < dp_cut = sqrt(dgn_ait e^{1.5 ln^2 sg} dgn_acc e^{1.5 ln^2 sg}) = 0.0745 um
           (modal_aero_rename)
  accum    otherwise
Per mode: N, dry volume, volume-weighted kappa, fixed sigma_g; number-median diameter from
N and V. Variant "mam4" applies the calcsize bounds (N adjusted so dgn stays within
[dgnumlo, dgnumhi]); "mam4_free" does not.

Usage: python3 tools/wp2_mam4_projection.py [runs/wp2] [--every 1]
Writes runs/wp2/mam4_projection.json and runs/wp2/mam4_projection.md.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.special import erfc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scrambox import case as casemod                       # noqa: E402
from scrambox.compile import build_grid                    # noqa: E402
from scrambox.partmc import _pmc_name, read_particles      # noqa: E402
from scrambox.post.ccn import (CPAIR, GRAVIT, LATVAP, MWH2O, PREF, R_UNIV, RAIR,  # noqa: E402
                               RH2O, RHOH2O, SURFTEN, TMELT, qsat_water)
from scrambox.post.common import KAPPA                     # noqa: E402

# CESM2 MAM4 physprops (mam4_mode1..4 files of the production cases): sigma_g, dgnumlo, dgnumhi (um)
MODES = {"accum": (1.6, 0.0535, 0.48), "aitken": (1.6, 0.0087, 0.052),
         "coarse": (1.2, 0.4, 40.0), "pcarbon": (1.6, 0.010, 0.100)}
DGN_NOMINAL = {"accum": 0.11, "aitken": 0.026}
DP_CUT = math.sqrt(DGN_NOMINAL["aitken"] * math.exp(1.5 * math.log(1.6) ** 2)
                   * DGN_NOMINAL["accum"] * math.exp(1.5 * math.log(1.6) ** 2))
D_COARSE_UM = 1.0
DR_AGE_UM = 8.0 * 4.76e-4
SUPERSATS = (0.1, 0.2, 0.5)
UPDRAFTS = (0.2, 0.5, 1.0)
SCENARIOS = ("urban_plume_aging", "seasalt_nitrate", "dust_crustal")
CARBON, MINERAL = ("BC", "POA"), ("OIN", "Na", "Cl")


def aten_of(T):
    return 2.0 * MWH2O * SURFTEN / (R_UNIV * T * RHOH2O)


def arg_thermo(T, P, w):
    """alpha*w, G and the eta prefactor of ndrop activate_modal (eta_m = pref / N_m)."""
    rho_air = P / (RAIR * T)
    diff0 = 0.211e-4 * (PREF / P) * (T / TMELT) ** 1.94
    conduct0 = (5.69 + 0.017 * (T - TMELT)) * 4.186e2 * 1e-5
    _, qs = qsat_water(T, P)
    dqsdt = LATVAP / (RH2O * T * T) * qs
    alpha = GRAVIT * (LATVAP / (CPAIR * RH2O * T * T) - 1.0 / (RAIR * T))
    gamma = (1.0 + LATVAP / CPAIR * dqsdt) / (rho_air * qs)
    G = 1.0 / (RHOH2O / (diff0 * rho_air * qs) + LATVAP * RHOH2O / (conduct0 * T) * (LATVAP / (RH2O * T) - 1.0))
    beta = 2.0 * math.pi * RHOH2O * G * gamma
    alw = alpha * w
    zeta = 2.0 / 3.0 * math.sqrt(alw) * aten_of(T) / math.sqrt(G)
    return zeta, alw ** 1.5 / (beta * math.sqrt(G))


def smax_arg(n, sc, f1, f2, zeta, pref):
    """ndrop maxsat for populations with number n (m-3), critical supersaturation sc."""
    ok = (n > 1e-20) & np.isfinite(sc)
    eta = pref / n[ok]
    g1 = np.sqrt(zeta / eta) ** 3
    g2 = np.sqrt(sc[ok] / np.sqrt(eta + 3.0 * zeta)) ** 3
    tot = np.sum((np.broadcast_to(f1, n.shape)[ok] * g1 + np.broadcast_to(f2, n.shape)[ok] * g2) / sc[ok] ** 2)
    return 1.0 / math.sqrt(tot) if tot > 0 else None


def particles(case_dir):
    """Per output: T, P, per-particle dry volumes (m3) by PartMC species, num conc (m-3)."""
    c = casemod.load(Path(case_dir) / "case.resolved.json")
    table = build_grid(c)["table"]
    dens = {_pmc_name(s): table.density_gcm3(s) * 1e3 for s in table.aerosol}
    kap = {_pmc_name(s): KAPPA.get(s, 0.1) for s in table.aerosol}
    for t, env, spec, mass, nconc, *_ in read_particles(case_dir):
        names = [s for s in spec if s != "H2O"]
        vol = np.array([mass[spec.index(s)] / dens[s] for s in names])
        ok = vol.sum(axis=0) > 0
        yield t, env.get("temperature", 298.0), env.get("pressure", 1e5), names, vol[:, ok], nconc[ok], kap


def truth_metrics(T, P, names, vol, n, kap):
    vt = vol.sum(axis=0)
    k = (vol * np.array([kap[s] for s in names])[:, None]).sum(axis=0) / vt
    r = (3 * vt / (4 * math.pi)) ** (1 / 3)
    A = aten_of(T)
    with np.errstate(divide="ignore", invalid="ignore"):
        sc = np.where(k > 0, 2 * A * np.sqrt(A / (27 * k * r ** 3)), np.inf)
    out = {"ccn": {str(s): float(n[sc <= s / 100].sum() / 1e6) for s in SUPERSATS}, "nd": {}, "smax": {}}
    for w in UPDRAFTS:
        zeta, pref = arg_thermo(T, P, w)
        smax = smax_arg(n, sc, 0.5, 1.0, zeta, pref)       # monodisperse populations (f1 0.5, f2 1)
        out["nd"][str(w)] = float(n[sc <= smax].sum() / 1e6) if smax else 0.0
        out["smax"][str(w)] = smax
    return out


def assign_modes(names, vol, kap):
    ix = {s: i for i, s in enumerate(names)}
    vt = vol.sum(axis=0)
    d_um = (6 * vt / math.pi) ** (1 / 3) * 1e6
    core = sum(vol[ix[s]] for s in CARBON if s in ix)
    mineral = sum(vol[ix[s]] for s in MINERAL if s in ix)
    shell = sum(vol[ix[s]] * kap[s] / KAPPA["PSO4"] for s in names if s not in CARBON + MINERAL)
    d_core_um = (6 * np.maximum(core, 0) / math.pi) ** (1 / 3) * 1e6
    need = math.pi * (d_core_um * 1e-6) ** 2 * DR_AGE_UM * 1e-6        # m3
    mode = np.where(d_um >= D_COARSE_UM, "coarse",
                    np.where((core > 0) & (mineral <= 0) & (shell < need), "pcarbon",
                             np.where(d_um < DP_CUT, "aitken", "accum")))
    return mode


def mam4_metrics(T, P, names, vol, n, kap, mode, bounded):
    A = aten_of(T)
    kv = np.array([kap[s] for s in names])
    pops = []
    for m, (sg, lo, hi) in MODES.items():
        sel = mode == m
        N = float(n[sel].sum())
        if N <= 0:
            continue
        V = float((vol[:, sel] * n[sel]).sum())
        kappa = float((vol[:, sel] * n[sel] * kv[:, None]).sum() / V)
        e45 = math.exp(4.5 * math.log(sg) ** 2)
        dgn = (6 * V / (math.pi * N * e45)) ** (1 / 3) * 1e6
        if bounded and not lo <= dgn <= hi:                 # calcsize: adjust number to the bound
            dgn = min(max(dgn, lo), hi)
            N = 6 * V / (math.pi * (dgn * 1e-6) ** 3 * e45)
        r3 = 3 * V / (4 * math.pi * e45 * N)                 # amcube (m3)
        sc = 2 * A * math.sqrt(A / (27 * max(kappa, 1e-30) * r3))
        pops.append({"mode": m, "N_cm3": N / 1e6, "dgn_um": dgn, "kappa": kappa, "sigma": sg, "sc": sc})
    N = np.array([p["N_cm3"] * 1e6 for p in pops])
    sc = np.array([p["sc"] for p in pops])
    lsg = np.log(np.array([p["sigma"] for p in pops]))
    frac = lambda s: 0.5 * erfc(2.0 / 3.0 * np.log(sc / s) / (math.sqrt(2) * lsg))
    out = {"modes": pops, "ccn": {str(s): float((N * frac(s / 100)).sum() / 1e6) for s in SUPERSATS},
           "nd": {}, "smax": {}}
    f1, f2 = 0.5 * np.exp(2.5 * lsg ** 2), 1.0 + 0.25 * lsg
    for w in UPDRAFTS:
        zeta, pref = arg_thermo(T, P, w)
        smax = smax_arg(N, sc, f1, f2, zeta, pref)
        out["nd"][str(w)] = float((N * frac(smax)).sum() / 1e6) if smax else 0.0
        out["smax"][str(w)] = smax
    return out


def run_scenario(sdir: Path, every: int):
    res = {}
    for pdir in sorted(sdir.glob("partmc_s[0-9]")):
        series = {"time_h": [], "truth": [], "mam4": [], "mam4_free": [], "mode_number_frac": []}
        for k, (t, T, P, names, vol, n, kap) in enumerate(particles(pdir)):
            if k % every:
                continue
            mode = assign_modes(names, vol, kap)
            series["time_h"].append(t / 3600)
            series["truth"].append(truth_metrics(T, P, names, vol, n, kap))
            series["mam4"].append(mam4_metrics(T, P, names, vol, n, kap, mode, True))
            series["mam4_free"].append(mam4_metrics(T, P, names, vol, n, kap, mode, False))
            series["mode_number_frac"].append({m: float(n[mode == m].sum() / n.sum()) for m in MODES})
        res[pdir.name] = series
        print(f"  {sdir.name}/{pdir.name}: {len(series['time_h'])} outputs", flush=True)
    return res


def scram_end(sdir: Path):
    """End-of-run CCN and Nd of the SCRAM runs and the PartMC->NC projections (summary.json)."""
    R = json.loads((sdir / "summary.json").read_text())["runs"]
    out = {}
    for key, r in R.items():
        if not isinstance(r, dict) or "ccn" not in r:
            continue
        out[key] = {"ccn": {s: r["ccn"][s][-1] for s in r["ccn"] if s in {str(x) for x in SUPERSATS}},
                    "nd": {w: r["nd"][w][-1] for w in r.get("nd", {}) if w in {str(x) for x in UPDRAFTS}}}
    return out


def mean_sd(vals):
    v = np.array([x for x in vals if x is not None], float)
    return (float(v.mean()), float(v.std())) if len(v) else (None, None)


def report(root: Path, allres: dict):
    L = ["# MAM4 vs SCRAM representation against PartMC (WP2)", "",
         "MAM4 = PartMC particles mapped onto the CESM2 MAM4 modes (CAM rules, calcsize bounds) and evaluated "
         "with CAM ccncalc / activate_modal; SCRAM NC1/NC6 = box runs; PartMC→NC = same particles on SCRAM "
         "bins and classes. Truth = particle-resolved. End of run, mean over seeds; bias vs truth in brackets.", ""]
    for scn, d in allres.items():
        seeds = list(d["mam4"].keys())
        tr = lambda q, x: mean_sd([d["mam4"][s]["truth"][-1][q][x] for s in seeds])[0]
        L += [f"## {scn}", "", "| quantity | truth | MAM4 | MAM4 (no size bounds) | SCRAM NC1 | SCRAM NC6 | PartMC→NC6 |",
              "|---|---|---|---|---|---|---|"]
        for q, xs, lab in (("ccn", SUPERSATS, "CCN S={}% (cm⁻³)"), ("nd", UPDRAFTS, "Nd w={} m/s (cm⁻³)")):
            for x in xs:
                t = tr(q, str(x))
                cells = []
                for v in (mean_sd([d["mam4"][s]["mam4"][-1][q][str(x)] for s in seeds])[0],
                          mean_sd([d["mam4"][s]["mam4_free"][-1][q][str(x)] for s in seeds])[0],
                          d["scram"].get("scram_NC1", {}).get(q, {}).get(str(x)),
                          d["scram"].get("scram_NC6", {}).get(q, {}).get(str(x)),
                          mean_sd([d["scram"].get(f"{s}_NC6", {}).get(q, {}).get(str(x)) for s in seeds])[0]):
                    cells.append("–" if v is None else f"{v:.4g} ({100 * (v / t - 1):+.0f}%)" if t else f"{v:.4g}")
                L.append(f"| {lab.format(x)} | {t:.4g} | " + " | ".join(cells) + " |")
        s0 = seeds[0]
        L += ["", "MAM4 modes at the end (seed 1): " + "; ".join(
            f"{p['mode']} N={p['N_cm3']:.3g} cm⁻³ dgn={p['dgn_um']:.3g} µm κ={p['kappa']:.3g}"
            for p in d["mam4"][s0]["mam4"][-1]["modes"]),
            "Number fraction by mode: " + ", ".join(f"{m} {v:.2f}" for m, v in d["mam4"][s0]["mode_number_frac"][-1].items()), ""]
    (root / "mam4_projection.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root", nargs="?", default="runs/wp2")
    ap.add_argument("--every", type=int, default=1)
    a = ap.parse_args()
    root = Path(a.root)
    allres = {}
    for scn in SCENARIOS:
        sdir = root / scn
        if not sdir.is_dir():
            continue
        allres[scn] = {"mam4": run_scenario(sdir, a.every), "scram": scram_end(sdir)}
    (root / "mam4_projection.json").write_text(json.dumps(allres))
    report(root, allres)


if __name__ == "__main__":
    main()
