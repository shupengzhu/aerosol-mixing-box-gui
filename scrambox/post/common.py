"""Shared pieces of the post-processors: the CESM-SCRAM species property tables
(kappa, optical type) and the per-(bin, class) population state of a box run.

Sources (CESM-SCRAM CAM component, carma_model_mod.F90):
  scram_species_kappa_value   kappa per species (MAM4 hygroscopicities)
  scram_optics_init           species -> optical type, refractive indices per type
  setupbins (explicit bins)   bin radius = geometric mean of the bin bounds
  scram_mam4_hygro_mod        MAM4 Koehler wet radius with hysteresis (rev3)
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

KAPPA = {"PMD": 0.068, "PBC": 1.0e-10, "PNA": 1.160, "PSO4": 0.507, "PNH4": 0.670,
         "PNO3": 0.670, "PHCL": 0.650, "PH2O": 0.0,
         "PSOA1P": 0.140, "PSOA2P": 0.140, "PSOA3P": 0.140, "PSOA4P": 0.140, "PSOA5P": 0.140,
         "PPOAlP": 1.0e-10}
IOPT_BC, IOPT_DUST, IOPT_SULF, IOPT_ORG, IOPT_SALT, IOPT_WATER = 1, 2, 3, 4, 5, 6
OPTYPE = {"PBC": IOPT_BC, "PMD": IOPT_DUST, "PSO4": IOPT_SULF, "PNH4": IOPT_SULF,
          "PNO3": IOPT_SULF, "PHCL": IOPT_SULF, "PNA": IOPT_SALT, "PH2O": IOPT_WATER}


def kappa_of(name: str, custom: dict | None = None) -> float:
    if custom and name in custom:
        return float(custom[name])
    if name in KAPPA:
        return KAPPA[name]
    raise KeyError(f"no kappa for species {name} (give one in the custom species definition)")


def optype_of(name: str) -> int:
    return OPTYPE.get(name, IOPT_ORG)          # SOA, POA and custom species -> organics


# ---------------------------------------------------------------- MAM4 water uptake
MAM4_RH_MAX, MAM4_RHCRYSTAL, MAM4_RHDELIQUES = 0.98, 0.35, 0.80


def mam4_kohler_rwet_m(rdry_m: float, hygro: float, s: float) -> float:
    eps, mw, surften, tair, rhow, ugascon = 1e-4, 18.0, 76.0, 273.0, 1.0, 8.3e7
    a = 2.0e4 * mw * surften / (ugascon * tair * rhow)          # um
    rdry = rdry_m * 1e6
    vol = rdry ** 3
    b = vol * max(hygro, 0.0)
    if vol <= 1e-12:
        return rdry_m
    ss = max(min(s, 1.0 - eps), 1e-10)
    slog = math.log(ss)
    p = (b / a) / (rdry * rdry)
    if p < eps:
        r = rdry * (1.0 + p / 3.0 / (1.0 - slog * rdry / a))
    else:
        g = lambda x: slog * x ** 4 - a * x ** 3 + (b - slog * vol) * x + a * vol
        rlo = rhi = rdry
        for _ in range(400):
            rhi *= 1.02
            if g(rhi) < 0.0:
                break
            rlo = rhi
        if g(rhi) >= 0.0:
            r = rdry
        else:
            for _ in range(60):
                rmid = 0.5 * (rlo + rhi)
                if g(rmid) > 0.0:
                    rlo = rmid
                else:
                    rhi = rmid
            r = 0.5 * (rlo + rhi)
    r = min(r, 30.0)
    return max(r, rdry) * 1e-6


def mam4_hysteresis_factor(rh: float) -> float:
    if rh < MAM4_RHCRYSTAL:
        return 0.0
    if rh < MAM4_RHDELIQUES:
        return (rh - MAM4_RHCRYSTAL) / max(1e-5, MAM4_RHDELIQUES - MAM4_RHCRYSTAL)
    return 1.0


def mam4_wet_radius_m(rdry_m: float, hygro: float, rh: float) -> float:
    r3dry = rdry_m ** 3
    r3wet = max(mam4_kohler_rwet_m(rdry_m, hygro, rh), rdry_m) ** 3
    wtr = max(r3wet - r3dry, 0.0) * mam4_hysteresis_factor(rh)
    return (r3dry + wtr) ** (1.0 / 3.0)


def kappa_wet_radius_cm(rdry_cm: float, rh: float, kappa_eff: float) -> float:
    """scram_kappa_wet_radius: MAM4 Koehler (273 K, RH <= 0.98, hysteresis)."""
    if rdry_cm <= 0:
        return rdry_cm
    rw = 1e2 * mam4_wet_radius_m(rdry_cm * 1e-2, max(kappa_eff, 0.0), min(max(rh, 0.0), MAM4_RH_MAX))
    return max(rw, rdry_cm)


# ---------------------------------------------------------------- population state
@dataclass
class Populations:
    """One output time of a run, as (bin, class) populations.
    number (nb, nc) #/m3; mass (nb, nc, naer) ug/m3; species/densities in kernel order."""
    number: np.ndarray
    mass: np.ndarray
    species: list
    density_gcm3: np.ndarray
    bounds_um: np.ndarray
    class_labels: list
    T: float
    RH: float
    P: float

    @property
    def nb(self):
        return self.number.shape[0]

    @property
    def nc(self):
        return self.number.shape[1]

    @property
    def iw(self):
        return self.species.index("PH2O")

    def r_mid_cm(self):
        """CARMA explicit-bin radius: geometric mean of the bin bounds (cm)."""
        return 0.5e-4 * np.sqrt(self.bounds_um[:-1] * self.bounds_um[1:])

    def r_low_up_cm(self):
        return 0.5e-4 * self.bounds_um[:-1], 0.5e-4 * self.bounds_um[1:]

    def internal_mean(self) -> "Populations":
        """All classes of a size bin merged into one internally mixed population."""
        return Populations(self.number.sum(axis=1, keepdims=True), self.mass.sum(axis=1, keepdims=True),
                           self.species, self.density_gcm3, self.bounds_um, ["internal mean"],
                           self.T, self.RH, self.P)


def populations(run: dict, it: int = -1) -> Populations:
    """Build the state at output index it from scrambox.run.load_run() / run_case() output."""
    g, r = run["grid"], run["results"]
    sp = g["species"]
    dens = np.array([float(sp["density_gcm3"][s]) for s in sp["aerosol"]])
    met = r["met"][it]
    return Populations(np.asarray(r["number_m3"][it], float), np.asarray(r["mass_ugm3"][it], float),
                       list(sp["aerosol"]), dens, np.asarray(g["bounds_um"], float),
                       list(g.get("class_labels", [])), float(met[0]), float(met[1]), float(met[2]))
