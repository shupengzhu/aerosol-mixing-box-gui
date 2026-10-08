"""Heterogeneous uptake (N2O5, HO2, NO3, NO2, ...) on a box state.

Rate per population p (CAM mo_usrrxt hetrxtrate, Fuchs-Sutugin-type resistance sum):
    k = sum_p  sfc_p / ( 0.5 dm_p / Dg + 4 / (c_gas gamma_p) )          [1/s]
with sfc in cm2/cm3, dm the wet diameter (cm), Dg = 0.1 cm2/s, c_gas = c0 sqrt(T) cm/s.

The surface follows scram_hetchem_sad_export (CESM-SCRAM hetfix): dry radius from the
population's number and dry mass, wet diameter from the population's own PH2O,
Dwet = Ddry (1 + V_H2O/V_dry)^(1/3).

gamma rules (case.postprocess.hetchem.rule, per gas):
  cesm        CESM: fixed gamma on the "chemically active" surface fraction
              f_chm = (PSO4 + PNH4 + PNO3 + PBC + PSOA1-5P) / dry mass  (exactly CESM)
  species     gamma_p = sum_i phi_i gamma_i over the dry volume fractions phi_i, with an
              editable per-species gamma table; full wet surface
  bertram_thornton (N2O5 only) Bertram & Thornton (2009) aqueous gamma from the
              population's water, nitrate and chloride molarities, optionally limited by an
              organic coating (Anttila et al., 2006; Riemer et al., 2009)
Each rule is applied to the class-resolved populations and to the internally mixed
(per size bin) reference, so the mixing-state effect on k and gamma_eff is explicit.
"""
from __future__ import annotations

import copy
import math

import numpy as np

from .common import Populations

DG_CM2_S = 0.1                                   # Dentener (1993), as CAM
C0 = {"N2O5": 1.40e3, "NO3": 1.85e3, "NO2": 2.15e3, "HO2": 2.53e3}   # c = c0 sqrt(T) cm/s (CAM)
GAMMA_CESM = {"N2O5": 0.10, "HO2": 0.20, "NO2": 1.0e-4, "NO3": 1.0e-3}
CHM_SPECIES = {"PSO4", "PNH4", "PNO3", "PBC", "PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P"}

# default per-species gamma table for the "species" rule; starts from the CESM values on
# the CESM chemically active species and zero elsewhere, i.e. reproduces CESM's surface
# selection. Users edit this table (rows species, columns gases) in the interface.
def default_species_gamma(species: list[str]) -> dict:
    tab = {}
    for g, v in GAMMA_CESM.items():
        tab[g] = {s: (v if s in CHM_SPECIES else 0.0) for s in species if s != "PH2O"}
    return tab


DEFAULTS = {"gases": ["N2O5", "HO2", "NO3", "NO2"],
            "rules": ["cesm", "species", "bertram_thornton"],
            "species_gamma": None,                  # None -> default_species_gamma
            "organic_coating": True,
            "bt_params": {"A_s": 3.2e-8, "beta_s1": 1.15e6, "delta_M1": 0.13,
                          "k3_k2b": 0.06, "k4_k2b": 29.0},
            "coating_params": {"H_aq_M_atm": 51.0, "D_aq_m2_s": 1.0e-9, "epsilon": 0.03}}

MW = {"PNO3": 62.0, "PHCL": 35.45, "PH2O": 18.015}
ORGANICS = {"PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P", "PPOAlP"}


def _geometry(pop: Populations):
    """Per population: number (#/cm3), dry and wet radius (cm), volumes (cm3/cm3)."""
    iw = pop.iw
    out = []
    for b in range(pop.nb):
        for c in range(pop.nc):
            n = pop.number[b, c]
            m = np.maximum(pop.mass[b, c], 0.0)
            mdry = m.sum() - m[iw]
            if n <= 0 or mdry <= 0:
                continue
            v = m / pop.density_gcm3 * 1e-12                     # (ug/m3)/(g/cm3) -> m3/m3
            vdry = v.sum() - v[iw]
            vwat = v[iw]
            rdry_m = (3.0 * vdry / (4.0 * math.pi * n)) ** (1.0 / 3.0)
            rwet_m = rdry_m * (1.0 + vwat / vdry) ** (1.0 / 3.0) if vwat > 0 else rdry_m
            out.append({"bin": b, "cls": c, "n_m3": n, "m": m, "mdry": mdry, "v": v,
                        "vdry": vdry, "vwat": vwat, "rdry_m": rdry_m, "rwet_m": rwet_m})
    return out


def _bt_gamma(p, species, prm, coat, cprm, c_ms, T):
    """Bertram & Thornton (2009) with an optional organic coating limit."""
    iw = species.index("PH2O")
    vwat_L = p["vwat"] * 1e3                                    # L per m3 air
    if vwat_L <= 0:
        return 0.0
    vsol_L = (p["vdry"] - _vol(p, species, ORGANICS | {"PBC", "PMD"})) * 1e3 + vwat_L
    h2o = p["m"][iw] * 1e-6 / MW["PH2O"] / vsol_L               # mol/L
    no3 = (p["m"][species.index("PNO3")] * 1e-6 / MW["PNO3"] / vsol_L) if "PNO3" in species else 0.0
    cl = (p["m"][species.index("PHCL")] * 1e-6 / MW["PHCL"] / vsol_L) if "PHCL" in species else 0.0
    k2f = prm["beta_s1"] - prm["beta_s1"] * math.exp(-prm["delta_M1"] * h2o)
    if no3 > 1e-12:
        frac = 1.0 - 1.0 / (prm["k3_k2b"] * h2o / no3 + 1.0 + prm["k4_k2b"] * cl / no3)
    else:
        frac = 1.0
    g_core = prm["A_s"] * k2f * frac
    if not coat:
        return g_core
    v_org = _vol(p, species, ORGANICS)
    if v_org <= 0 or g_core <= 0:
        return g_core
    rp = p["rwet_m"]
    rc = rp * max(1.0 - v_org / (p["vdry"] + p["vwat"]), 0.0) ** (1.0 / 3.0)
    l = max(rp - rc, 1e-12)
    R_atm = 0.082057                                            # L atm / (mol K)
    hd = cprm["epsilon"] * cprm["H_aq_M_atm"] * cprm["D_aq_m2_s"]   # M/atm m2/s
    g_coat = 4.0 * R_atm * T * hd * rc / (c_ms * l * rp)
    return 1.0 / (1.0 / g_core + 1.0 / g_coat)


def _vol(p, species, names):
    return sum(p["v"][i] for i, s in enumerate(species) if s in names)


def _rate(pops, gamma_fn, area_fn, c_cms):
    """k (1/s) and surface (cm2/cm3) summed over populations; gamma via callback."""
    k = sfc_tot = gsfc = 0.0
    for p in pops:
        sfc = area_fn(p)                                        # cm2/cm3
        if sfc <= 0:
            continue
        g = gamma_fn(p)
        dm_cm = 2.0 * p["rwet_m"] * 1e2
        sfc_tot += sfc
        gsfc += g * sfc
        if g > 0:
            k += sfc / (0.5 * dm_cm / DG_CM2_S + 4.0 / (c_cms * g))
    return k, sfc_tot, gsfc


def compute(pop: Populations, opts: dict | None = None) -> dict:
    o = copy.deepcopy(DEFAULTS)
    o.update(opts or {})
    sp = pop.species
    table = o["species_gamma"] or default_species_gamma(sp)
    T = pop.T
    out = {"T": T, "rules": {}, "sad_um2_cm3": None}
    for label, state in (("resolved", pop), ("internal_mean", pop.internal_mean())):
        pops = _geometry(state)
        wet_area = lambda p: 1e-2 * p["n_m3"] * math.pi * (2 * p["rwet_m"]) ** 2     # cm2/cm3
        chm_area = lambda p: wet_area(p) * sum(p["m"][i] for i, s in enumerate(sp) if s in CHM_SPECIES) / p["mdry"]
        if label == "resolved":
            out["sad_um2_cm3"] = float(sum(wet_area(p) for p in pops) * 1e8)
            out["sad_by_class_um2_cm3"] = [float(sum(wet_area(p) for p in pops if p["cls"] == c) * 1e8)
                                          for c in range(state.nc)]
        for gas in o["gases"]:
            c_cms = C0[gas] * math.sqrt(T)
            res = {}
            for rule in o["rules"]:
                if rule == "cesm":
                    k, s, gs = _rate(pops, lambda p: GAMMA_CESM[gas], chm_area, c_cms)
                elif rule == "species":
                    gt = table.get(gas, {})

                    def gfun(p):
                        vd = p["vdry"]
                        return sum(p["v"][i] / vd * float(gt.get(s_, 0.0)) for i, s_ in enumerate(sp) if s_ != "PH2O")
                    k, s, gs = _rate(pops, gfun, wet_area, c_cms)
                elif rule == "bertram_thornton":
                    if gas != "N2O5":
                        continue
                    k, s, gs = _rate(pops, lambda p: _bt_gamma(p, sp, o["bt_params"], o["organic_coating"],
                                                              o["coating_params"], c_cms * 1e-2, T),
                                     wet_area, c_cms)
                else:
                    raise ValueError(f"unknown hetchem rule {rule}")
                # effective gamma: the gamma that gives k on the total wet surface (free molecular)
                sad_all = sum(wet_area(p) for p in pops)
                g_eff = 4.0 * k / (c_cms * sad_all) if sad_all > 0 else 0.0
                res[rule] = {"k_s1": k, "lifetime_h": (1.0 / k / 3600.0 if k > 0 else None),
                             "gamma_eff": g_eff, "gamma_surface_mean": (gs / s if s > 0 else 0.0)}
            out["rules"].setdefault(gas, {})[label] = res
    out["options"] = {k: v for k, v in o.items() if k != "species_gamma"}
    out["species_gamma"] = table
    return out


def series(run: dict, opts: dict | None = None, every: int = 1) -> dict:
    from .common import populations
    times = run["results"]["time_s"]
    idx = list(range(0, len(times), max(1, every)))
    if idx[-1] != len(times) - 1:
        idx.append(len(times) - 1)
    out = {"time_s": [float(times[i]) for i in idx], "states": []}
    for i in idx:
        out["states"].append(compute(populations(run, i), opts))
    return out
