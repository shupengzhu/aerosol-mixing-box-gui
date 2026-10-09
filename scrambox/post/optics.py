"""Aerosol optics of a box state, identical to CESM-SCRAM CARMA_DiagnoseBulk.

Options (case.postprocess.optics, same meaning as the CESM cfg keys):
  mixing               0 per-type Mie tables (external), 1 volume-mixed refractive index
                       (online Mie), 2 reduced-shell chi_w LUT (production),
                       3 coated core-shell LUT
  perclass             True: each composition class separately (production); False: classes
                       of a bin summed first
  prognostic_moments   False (production): dry radius = bin reference radius and particle
                       mass from the dry-mixture density; True: dry radius from N and M
  dryshell_selector    0 sulfate-like (production), 1 BrC-like, 2 organic fraction > 0.5
RH is the clear-sky ambient RH capped at 0.98 (Optics_water = 0, Optics_rwet = 0), and the
wet radius comes from the same MAM4 kappa-Koehler closure as CESM.

CESM accumulates col_mass * ext_per_mass; in a box the column mass is replaced by the dry
mass concentration, so tau becomes an extinction coefficient (1/m, reported in 1/Mm).

Known CESM behaviours reproduced on purpose (non-production paths):
  mixing 0  the RH node index is computed as if mie_rh were uniform (rh*(NMIE_RH-1)) while
            the grid is 0,.5,.7,.8,.9,.95,.98,.99, so RH 0.82 reads the 0.95-0.98 tables
            (growth strongly overestimated);
  mixing 1  the volume-mixed refractive index excludes water but is applied to the wet
            radius, so absorption by BC-containing particles is overestimated.
"""
from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

import numpy as np

from .common import (IOPT_BC, IOPT_DUST, IOPT_ORG, IOPT_SALT, IOPT_SULF, IOPT_WATER,
                     Populations, kappa_of, kappa_wet_radius_cm, optype_of)
from .mie import mie_q

LUT_DIR = Path(__file__).resolve().parents[2] / "kernel" / "data" / "lut"
RS_LUT = "scram_reducedshell_optics_lut_cdf5.nc"
CS_LUT = "scram_coreshell_optics_lut.nc"
NLW, NSW = 16, 14
IDX_SW_DIAG = 10                       # RRTMG SW band containing 550 nm (radconstants)
MIE_RH = np.array([0.0, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99])
NTYPE = 6

# Refractive index per optical type and CARMA wavelength (carma_model_mod scram_optics_init)
_RI = {
    IOPT_BC: ([2.50, 2.30, 2.20, 2.10, 2.05, 2.00, 1.98, 1.96, 1.93, 1.91, 1.90, 1.89, 1.88, 1.87, 1.86, 1.85] + [1.85] * 14,
              [1.20, 1.10, 1.05, 1.00, 0.95, 0.90, 0.87, 0.85, 0.82, 0.80, 0.78, 0.76, 0.75, 0.74, 0.73, 0.72] + [0.71] * 14),
    IOPT_DUST: ([2.40, 2.00, 2.20, 2.00, 1.90, 1.80, 1.60, 1.50, 1.80, 1.90, 2.00, 2.10, 2.00, 1.90, 1.90, 1.80,
                 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.55, 1.56, 1.57, 1.53],
                [0.60, 0.50, 0.50, 0.60, 0.50, 0.70, 0.90, 1.00, 0.30, 0.20, 0.15, 0.10, 0.20, 0.50, 0.30, 0.10,
                 0.008, 0.008, 0.008, 0.008, 0.008, 0.006, 0.005, 0.005, 0.004, 0.005, 0.007, 0.010, 0.013, 0.008]),
    IOPT_SULF: ([1.82, 1.72, 1.62, 1.48, 1.50, 1.96, 1.66, 1.50, 1.52, 1.42, 1.82, 1.59, 1.53, 1.47, 1.46, 1.45] + [1.52] * 14,
                [0.20, 0.15, 0.30, 0.25, 0.15, 0.40, 0.35, 0.05, 0.10, 0.30, 0.60, 0.03, 0.01, 0.01, 0.01, 0.005] + [1.0e-7] * 14),
    IOPT_ORG: ([1.80, 1.68, 1.60, 1.55, 1.54, 1.53, 1.52, 1.50, 1.49, 1.48, 1.47, 1.46, 1.45, 1.44, 1.44, 1.43,
                1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.53, 1.54, 1.54, 1.53],
               [0.20, 0.13, 0.20, 0.15, 0.10, 0.08, 0.06, 0.04, 0.03, 0.02, 0.02, 0.015, 0.01, 0.01, 0.008, 0.005,
                0.006, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.006, 0.008, 0.010, 0.015, 0.020, 0.006]),
    IOPT_SALT: ([1.73, 1.60, 1.54, 1.51, 1.50, 1.49, 1.48, 1.48, 1.50, 1.40, 1.30, 1.48, 1.50, 1.52, 1.52, 1.52] + [1.50] * 14,
                [0.25, 0.15, 0.02, 0.005, 0.005, 0.003, 0.002, 0.001, 0.005, 0.20, 0.10, 0.001, 0.005, 0.008, 0.003, 0.001] + [1.0e-7] * 14),
    IOPT_WATER: ([1.532, 1.524, 1.420, 1.274, 1.161, 1.142, 1.232, 1.266, 1.296, 1.321, 1.342, 1.315, 1.330, 1.339, 1.350, 1.408,
                  1.408, 1.324, 1.277, 1.302, 1.312, 1.321, 1.323, 1.327, 1.331, 1.334, 1.340, 1.349, 1.362, 1.291],
                 [0.336, 0.360, 0.426, 0.403, 0.321, 0.115, 0.047, 0.039, 0.034, 0.034, 0.092, 0.012, 0.013, 0.010, 0.005, 0.014,
                  0.014, 0.158, 1.52e-3, 1.16e-3, 2.36e-4, 1.71e-4, 2.44e-5, 3.13e-6, 3.74e-8, 1.64e-9, 2.49e-9, 1.52e-8, 3.35e-8, 0.034]),
}
REFIDX = np.zeros((30, NTYPE + 1), dtype=complex)
for _t, (_re, _im) in _RI.items():
    REFIDX[:, _t] = np.array(_re) + 1j * np.array(_im)

DEFAULTS = {"mixing": 2, "perclass": True, "prognostic_moments": False, "dryshell_selector": 0}


# ---------------------------------------------------------------- LUTs
@lru_cache(maxsize=2)
def _lut(name):
    import netCDF4
    path = LUT_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"{path} missing - run kernel/fetch_luts.sh")
    from ..ncio import NC_LOCK
    with NC_LOCK, netCDF4.Dataset(path) as d:
        out = {k: np.array(d[k][:], dtype=float) for k in ("wavelength", "r_wet", "eta", "fbc", "Qext", "Qsca", "asym")}
        if "chi_w" in d.variables:
            out["chi_w"] = np.array(d["chi_w"][:], dtype=float)
            if int(d.getncattr("failed_nodes")) != 0:
                raise ValueError("reduced-shell LUT contains failed nodes")
    out["log_r"] = np.log(out["r_wet"])
    return out


def wavelengths_um() -> np.ndarray:
    """CARMA / RRTMG band wavelengths (the reduced-shell LUT grid, validated in CESM)."""
    return _lut(RS_LUT)["wavelength"] * 1e4


def _bracket(arr, val):
    """find_bracket of the CESM optics modules (clamped)."""
    n = arr.size
    if val <= arr[0]:
        return 0, 0.0
    if val >= arr[-1]:
        return n - 2, 1.0
    i = int(np.searchsorted(arr, val, side="right") - 1)
    i = min(max(i, 0), n - 2)
    w = (val - arr[i]) / (arr[i + 1] - arr[i]) if arr[i + 1] > arr[i] else 0.0
    return i, w


def _interp_rs(rwet_cm, eta, fbc, chi, idry):
    """scram_reducedshell_interp for all 30 wavelengths: arrays (30,) Qext, Qsca, g."""
    L = _lut(RS_LUT)
    ir, wr = _bracket(L["log_r"], math.log(max(rwet_cm, L["r_wet"][0])))
    ie, we = _bracket(L["eta"], min(max(eta, 0.0), 1.0))
    ifb, wf = _bracket(L["fbc"], min(max(fbc, 0.0), 1.0))
    ic, wc = _bracket(L["chi_w"], min(max(chi, L["chi_w"][0]), L["chi_w"][-1]))
    res = []
    for key in ("Qext", "Qsca", "asym"):
        F = L[key]
        v = np.zeros(F.shape[0])
        for dc, wcl in ((0, 1 - wc), (1, wc)):
            for df, wfl in ((0, 1 - wf), (1, wf)):
                for de, wel in ((0, 1 - we), (1, we)):
                    for dr, wrl in ((0, 1 - wr), (1, wr)):
                        v = v + wcl * wfl * wel * wrl * F[:, ir + dr, ie + de, ifb + df, ic + dc, idry - 1]
        res.append(v)
    qe, qs, g = res
    qe = np.maximum(qe, 0.0)
    return qe, np.minimum(np.maximum(qs, 0.0), qe), np.clip(g, -1.0, 1.0)


def _interp_cs(rwet_cm, eta, fbc, ishell):
    """scram_coreshell_interp (nearest shell class, trilinear r/eta/fbc)."""
    L = _lut(CS_LUT)
    ish = min(max(ishell, 1), 5) - 1
    ir, wr = _bracket(L["log_r"], math.log(max(rwet_cm, L["r_wet"][0])))
    ie, we = _bracket(L["eta"], min(max(eta, 0.0), 1.0))
    ifb, wf = _bracket(L["fbc"], min(max(fbc, 0.0), 1.0))
    res = []
    for key in ("Qext", "Qsca", "asym"):
        F = L[key][:, :, :, :, ish]
        q00 = (1 - wf) * F[:, ir, ie, ifb] + wf * F[:, ir, ie, ifb + 1]
        q01 = (1 - wf) * F[:, ir, ie + 1, ifb] + wf * F[:, ir, ie + 1, ifb + 1]
        q10 = (1 - wf) * F[:, ir + 1, ie, ifb] + wf * F[:, ir + 1, ie, ifb + 1]
        q11 = (1 - wf) * F[:, ir + 1, ie + 1, ifb] + wf * F[:, ir + 1, ie + 1, ifb + 1]
        res.append((1 - wr) * ((1 - we) * q00 + we * q01) + wr * ((1 - we) * q10 + we * q11))
    qe, qs, g = res
    qe = np.maximum(qe, 0.0)
    return qe, np.minimum(np.maximum(qs, 0.0), qe), np.clip(g, -1.0, 1.0)


def _classify_shell(v_sulf, v_org, v_salt, v_water):
    """scram_coreshell_classify."""
    vs = v_sulf + v_org + v_salt + v_water
    if vs < 1e-30:
        return 4
    if v_water / vs > 0.60:
        return 5
    vd = v_sulf + v_org + v_salt
    if vd < 1e-30:
        return 5
    if v_salt / vd > 0.40:
        return 2
    if v_org / vd > 0.40:
        return 3
    if v_sulf / vd > 0.30:
        return 1
    return 4


@lru_cache(maxsize=8)
def _type_tables(bounds_key, rho_ref):
    """Optics_mixing 0: per-type Mie tables (bin, RH node, wavelength, type) in m2/kg,
    built like scram_optics_init (reference density = PBC density, the number carrier)."""
    bounds = np.array(bounds_key)
    r = 0.5e-4 * np.sqrt(bounds[:-1] * bounds[1:])
    rmass_g = 4.0 / 3.0 * math.pi * rho_ref * r ** 3
    wave = wavelengths_um() * 1e-4
    nb = r.size
    ext = np.zeros((nb, MIE_RH.size, 30, NTYPE + 1))
    ssa = np.zeros_like(ext)
    asm = np.zeros_like(ext)
    for it in range(1, NTYPE + 1):
        kap = {IOPT_BC: kappa_of("PBC"), IOPT_DUST: kappa_of("PMD"), IOPT_SULF: kappa_of("PSO4"),
               IOPT_ORG: kappa_of("PSOA1P"), IOPT_SALT: kappa_of("PNA")}.get(it, 0.0)
        rw = np.array([[kappa_wet_radius_cm(r[b], MIE_RH[k], kap) for k in range(MIE_RH.size)] for b in range(nb)])
        for iw in range(30):
            x = 2 * math.pi * rw.ravel() / wave[iw]
            qe, qs, g = mie_q(x, np.full(x.size, REFIDX[iw, it]))
            qe, qs, g = (a.reshape(rw.shape) for a in (qe, qs, g))
            area = math.pi * (rw * 1e-2) ** 2
            if iw < NLW:
                ext[:, :, iw, it] = (qe - qs) * area / (rmass_g[:, None] * 1e-3)
            else:
                ext[:, :, iw, it] = qe * area / (rmass_g[:, None] * 1e-3)
                ssa[:, :, iw, it] = np.where(qe > 0, qs / np.where(qe > 0, qe, 1), 0)
                asm[:, :, iw, it] = g
    return ext, ssa, asm


def compute(pop: Populations, opts: dict | None = None, custom_kappa: dict | None = None) -> dict:
    """Optical properties of one state. Returns per-band arrays (1/Mm) and diagnostics."""
    o = dict(DEFAULTS, **(opts or {}))
    mixing, perclass = int(o["mixing"]), bool(o["perclass"])
    wave = wavelengths_um()
    sp = pop.species
    kap = np.array([kappa_of(s, custom_kappa) for s in sp])
    otype = np.array([optype_of(s) for s in sp])
    rho = pop.density_gcm3
    iw_h2o = pop.iw
    rh = min(max(pop.RH, 0.0), 0.98)                        # clear-sky MAM4 RH, cloud fraction 0
    r_ref = pop.r_mid_cm()
    tau = np.zeros(30)                                      # ext (LW: absorption) [1/m]
    sig = np.zeros(30)
    asy = np.zeros(30)
    bc_ext = dust_ext = 0.0
    ext_cls = np.zeros((pop.nb, pop.nc if perclass else 1))
    abs_bc_coated = abs_bc_bare = 0.0
    if mixing == 0:
        tab_ext, tab_ssa, tab_asm = _type_tables(tuple(pop.bounds_um), float(rho[sp.index("PBC")]) if "PBC" in sp else 2.25)
        if rh <= MIE_RH[0]:
            irh, wrh = 0, 0.0
        else:
            irh = min(int(math.floor(rh * (MIE_RH.size - 1))), MIE_RH.size - 2)
            wrh = min(max(rh * (MIE_RH.size - 1) - irh, 0.0), 1.0)
    number = pop.number if perclass else pop.number.sum(axis=1, keepdims=True)
    mass = pop.mass if perclass else pop.mass.sum(axis=1, keepdims=True)
    for b in range(pop.nb):
        for c in range(number.shape[1]):
            m = np.maximum(mass[b, c], 0.0)                   # ug/m3
            dry = np.ones(len(sp), bool); dry[iw_h2o] = False
            mdry = m[dry].sum()
            if mdry <= 0.0:
                continue
            vol = np.where(dry, m / rho, 0.0)                 # proportional to volume
            vtot = vol.sum()
            mass_kgm3 = mdry * 1e-9
            # particle mass and dry radius (scram_optics_reference_mass / _moments)
            if o["prognostic_moments"]:
                if number[b, c] <= 0:
                    continue
                pm_kg = mdry * 1e-9 / number[b, c]
                vcm3 = pm_kg * 1e3 * (vtot / mdry)
                rdry = (vcm3 * 3.0 / (4.0 * math.pi)) ** (1.0 / 3.0)
            else:
                rdry = r_ref[b]
                pm_kg = (4.0 * math.pi / 3.0) * rdry ** 3 / (vtot / mdry) * 1e-3
            keff = float((kap[dry] * vol[dry]).sum() / vtot)
            rwet = kappa_wet_radius_cm(rdry, rh, keff)
            area = math.pi * (rwet * 1e-2) ** 2
            v_bc = vol[otype == IOPT_BC].sum()
            v_du = vol[otype == IOPT_DUST].sum()
            if mixing == 0:
                for j, s in enumerate(sp):
                    if j == iw_h2o or m[j] <= 0:
                        continue
                    t = otype[j]
                    e = (1 - wrh) * tab_ext[b, irh, :, t] + wrh * tab_ext[b, irh + 1, :, t]
                    s_ = (1 - wrh) * tab_ssa[b, irh, :, t] + wrh * tab_ssa[b, irh + 1, :, t]
                    g_ = (1 - wrh) * tab_asm[b, irh, :, t] + wrh * tab_asm[b, irh + 1, :, t]
                    cm = m[j] * 1e-9
                    tau += cm * e
                    sig[NLW:] += (cm * e * s_)[NLW:]
                    asy[NLW:] += (cm * e * s_ * g_)[NLW:]
                    if t == IOPT_BC:
                        bc_ext += cm * e[NLW + IDX_SW_DIAG - 1]
                    if t == IOPT_DUST:
                        dust_ext += cm * e[NLW + IDX_SW_DIAG - 1]
                    ext_cls[b, c] += cm * e[NLW + IDX_SW_DIAG - 1]
                continue
            if mixing == 1:
                ri = (REFIDX[:, otype[dry]] * vol[dry]).sum(axis=1) / vtot
                qe, qs, g = mie_q(2 * math.pi * rwet / (wave * 1e-4), ri)
            elif mixing == 2:
                v_core = v_bc + v_du
                v_shell_dry = vtot - v_core
                v_org = vol[(otype == IOPT_ORG) & dry].sum()
                swell = max((rwet / rdry) ** 3, 1.0) if rdry > 0 else 1.0
                v_eff = vtot * swell
                v_shell_wet = max(v_eff - v_core, v_shell_dry)
                eta = (v_core / v_eff) ** (1 / 3) if v_core > 1e-30 * v_eff else 0.0
                fbc = v_bc / v_core if v_core > 1e-30 else 0.0
                chi = (v_shell_wet - v_shell_dry) / v_shell_wet if v_shell_wet > 1e-30 else 0.0
                chi = min(max(chi, 0.0), 0.99)
                sel = int(o["dryshell_selector"])
                if sel == 0:
                    idry = 1
                elif sel == 1:
                    idry = 2
                else:
                    idry = 2 if (v_shell_dry > 1e-30 and v_org / v_shell_dry > 0.5) else 1
                qe, qs, g = _interp_rs(rwet, eta, fbc, chi, idry)
                v_frac_tot = v_eff
            else:   # mixing 3
                v_s = vol[otype == IOPT_SULF].sum()
                v_o = vol[(otype == IOPT_ORG) & dry].sum()
                v_sa = vol[otype == IOPT_SALT].sum()
                v_core = v_bc + v_du
                eta = (v_core / vtot) ** (1 / 3) if v_core > 1e-30 * vtot else 0.0
                fbc = v_bc / v_core if v_core > 1e-30 else 0.0
                qe, qs, g = _interp_cs(rwet, eta, fbc, _classify_shell(v_s, v_o, v_sa, 0.0))
                v_frac_tot = vtot
            e = qe * area / pm_kg                               # m2/kg
            ssa = np.where(qe > 0, qs / np.where(qe > 0, qe, 1), 0.0)
            tau[:NLW] += mass_kgm3 * np.maximum(qe - qs, 0.0)[:NLW] * area / pm_kg
            tau[NLW:] += mass_kgm3 * e[NLW:]
            sig[NLW:] += mass_kgm3 * (e * ssa)[NLW:]
            asy[NLW:] += mass_kgm3 * (e * ssa * g)[NLW:]
            i550 = NLW + IDX_SW_DIAG - 1
            vt = vtot if mixing == 1 else v_frac_tot
            bc_ext += mass_kgm3 * e[i550] * v_bc / max(vt, 1e-30)
            dust_ext += mass_kgm3 * e[i550] * v_du / max(vt, 1e-30)
            ext_cls[b, c] += mass_kgm3 * e[i550]
            if v_bc > 0:
                # absorption enhancement of BC: this population vs its BC as a bare sphere
                abs_bc_coated += mass_kgm3 * (qe[i550] - qs[i550]) * area / pm_kg
                rbc = rdry * (v_bc / vtot) ** (1 / 3)
                qeb, qsb, _ = mie_q([2 * math.pi * rbc / (wave[i550] * 1e-4)], [REFIDX[i550, IOPT_BC]])
                n_part = mass_kgm3 / pm_kg
                abs_bc_bare += n_part * (qeb[0] - qsb[0]) * math.pi * (rbc * 1e-2) ** 2
    sw = slice(NLW, 30)
    ext_sw = tau[sw]
    ssa_sw = np.where(ext_sw > 0, sig[sw] / np.where(ext_sw > 0, ext_sw, 1), 0.0)
    g_sw = np.where(sig[sw] > 0, asy[sw] / np.where(sig[sw] > 0, sig[sw], 1), 0.0)
    i = IDX_SW_DIAG - 1
    to_Mm = 1e6
    return {
        "wavelength_sw_um": wave[sw].tolist(), "wavelength_lw_um": wave[:NLW].tolist(),
        "ext_sw_Mm": (ext_sw * to_Mm).tolist(), "ssa_sw": ssa_sw.tolist(), "g_sw": g_sw.tolist(),
        "abs_lw_Mm": (tau[:NLW] * to_Mm).tolist(),
        "ext550_Mm": float(ext_sw[i] * to_Mm), "ssa550": float(ssa_sw[i]), "g550": float(g_sw[i]),
        "abs550_Mm": float(ext_sw[i] * (1 - ssa_sw[i]) * to_Mm),
        "bc_ext550_Mm": float(bc_ext * to_Mm), "dust_ext550_Mm": float(dust_ext * to_Mm),
        "ext550_by_class_Mm": (ext_cls.sum(axis=0) * to_Mm).tolist(),
        "ext550_by_bin_Mm": (ext_cls.sum(axis=1) * to_Mm).tolist(),
        "bc_abs_enhancement": (float(abs_bc_coated / abs_bc_bare) if abs_bc_bare > 0 else None),
        "options": o,
    }


def series(run: dict, opts: dict | None = None, every: int = 1) -> dict:
    """Optics at every `every`-th output time, plus the internally mixed reference."""
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
