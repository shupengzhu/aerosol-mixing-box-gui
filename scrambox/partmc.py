"""PartMC particle-resolved reference runs from the same case.json.

PartMC (with the pmc_scram_thermo module: SCRAM's dynamic ISORROPIA scheme per particle,
CESM VBS SOA bulk equilibrium, SOA photolysis) is driven by the identical scenario:
species, gases, initial and background states, emissions and their time profiles,
dilution, temperature, pressure, RH and J(NO2). PartMC's own water condensation and
nucleation are off; coagulation (Brownian) follows the case.

  write_partmc(case, run_dir)   -> PartMC input files + thermo namelist
  run_partmc(case, run_dir)     -> runs PartMC, then project() + truth()
  project(run_dir)              -> PartMC particles on the SCRAM bin x class grid of the
                                   case, saved like a SCRAM run (results.npz, grid.json, ...)
                                   so post-processing and comparisons apply unchanged
  truth(run_dir)                -> particle-resolved metrics (CCN spectrum from per-particle
                                   kappa-Koehler, mixing-state index chi, BC coating)

Units: PartMC height is fixed at 1 m, so its per-m2 emissions are per-m3 like the box.
"""
from __future__ import annotations

import glob
import json
import math
import os
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np

from . import case as casemod
from .compile import build_grid, met_table
from .grid import find_class
from .post.common import KAPPA

PARTMC_BIN = {"v23": os.environ.get("PARTMC_V23", "/data/shupengz/Model/partmc/build_scram_v23/partmc"),
              "v1": os.environ.get("PARTMC_V1", "/data/shupengz/Model/partmc/build_scram_v1scram/partmc")}

AERO_MAP = {"PMD": "OIN", "PBC": "BC", "PNA": "Na", "PSO4": "SO4", "PNH4": "NH4", "PNO3": "NO3",
            "PHCL": "Cl", "PSOA1P": "SOA1", "PSOA2P": "SOA2", "PSOA3P": "SOA3", "PSOA4P": "SOA4",
            "PSOA5P": "SOA5", "PPOAlP": "POA", "PH2O": "H2O"}
GAS_MAP = {"SULF": "H2SO4", "NH3": "NH3", "HNO3": "HNO3", "HCL": "HCl", "SOA1P": "SOAG0",
           "SOA2P": "SOAG1", "SOA3P": "SOAG2", "SOA4P": "SOAG3", "SOA5P": "SOAG4"}
GAS_MW = {"SULF": 98.0, "NH3": 17.0, "HNO3": 63.0, "HCL": 36.5, "SOA1P": 250.0, "SOA2P": 250.0,
          "SOA3P": 250.0, "SOA4P": 250.0, "SOA5P": 250.0}
R = 8.314472


def _pmc_name(s):
    return AERO_MAP.get(s, s)


def _ppb(ug, mw, T, P):
    return ug / (mw * 1e6) / (P / (R * T)) * 1e9


def _series_file(path, t, v):
    Path(path).write_text("".join(f"{a:.6f} {b:.8e}\n" for a, b in zip(t, v)))


def _modes(spec, table, key_n, prefix, d):
    """case particle spec -> PartMC aero_dist text (+ composition files)."""
    lines = []
    if spec.get("mode", "lognormal") == "lognormal":
        for k, m in enumerate(spec.get("modes", [])):
            name = f"{prefix}{k}"
            fr = {_pmc_name(s): float(v) for s, v in (m.get("mass_fraction") or {}).items() if float(v) > 0}
            (d / f"{name}_comp.dat").write_text("".join(f"{s} {v}\n" for s, v in fr.items()))
            lines.append(f"mode_name {name}\nmass_frac {name}_comp.dat\ndiam_type geometric\nmode_type log_normal\n"
                         f"num_conc {float(m[key_n]) * 1e6:.6e}\ngeom_mean_diam {float(m['dg_um']) * 1e-6:.6e}\n"
                         f"log10_geom_std_dev {math.log10(float(m['sigma_g'])):.6f}\n")
    elif spec.get("mode") == "bin_mass":
        raise ValueError("PartMC reference runs need log-normal particle modes (bin_mass not supported yet)")
    elif spec.get("mode") == "table":
        raise ValueError("PartMC reference runs need log-normal particle modes (table not supported yet)")
    if not lines:
        return "# no particles\n"
    return "\n".join(lines)


def write_partmc(user_case: dict, run_dir, n_part: int = 10000, del_t: float = 60.0, seed: int = 1,
                 weighting_exponent: float = -1.0) -> dict:
    c = casemod.resolve(user_case)
    if c["kernel"]["isorropia"] == "none":
        raise ValueError("PartMC reference runs use SCRAM's ISORROPIA scheme: choose kernel v23 or v1")
    d = Path(run_dir)
    d.mkdir(parents=True, exist_ok=True)
    g = build_grid(c)
    table = g["table"]
    met = met_table(c)                      # rows: t, T, RH, P, J, efac, gfac, kdil, n2o5
    times = met[0]
    T0, P0 = float(met[1][0]), float(met[3][0])

    # species and gases
    rows = []
    for s in table.aerosol:
        rho = table.density_gcm3(s) * 1000.0
        mw = table.mw_gmol(s) / 1000.0
        kap = KAPPA.get(s, 0.1)
        rows.append(f"{_pmc_name(s):8s} {rho:8.1f} 0 {mw:.6e} {kap:.4g} 0 0")
    (d / "aero_data.dat").write_text("# name dens(kg/m3) ions MW(kg/mol) kappa abifm_m abifm_c\n" + "\n".join(rows) + "\n")
    gases = [s for s in GAS_MAP if s in table.gases]
    (d / "gas_data.dat").write_text("".join(GAS_MAP[s] + "\n" for s in gases))
    g0 = c["initial"].get("gas_ugm3") or {}
    (d / "gas_init.dat").write_text("# species ppb\n" + "".join(
        f"{GAS_MAP[s]} {_ppb(float(g0.get(s, 0.0)), GAS_MW[s], T0, P0):.6e}\n" for s in gases))

    # particles: initial, background (= initial unless given), emissions
    (d / "aero_init_dist.dat").write_text(_modes(c["initial"]["particles"], table, "number_cm3", "init", d))
    bg = (c["environment"].get("dilution") or {}).get("background", "initial")
    if bg in (None, "initial"):
        shutil.copy(d / "aero_init_dist.dat", d / "aero_back_dist.dat")
        gb = g0
    else:
        (d / "aero_back_dist.dat").write_text(_modes(bg.get("particles") or {"modes": []}, table, "number_cm3", "back", d))
        gb = bg.get("gas_ugm3") or {}
    (d / "aero_emit_dist.dat").write_text(_modes(c["emissions"]["particles"], table, "number_cm3_s", "emit", d))

    # PartMC holds emission/background rates constant from one time point to the next;
    # the box applies them per step with the value at the step start. Writing the rates
    # on the box step grid with those start values gives identical source integrals.
    dt_box = float(c["time"]["dt_s"])
    st = np.arange(0.0, float(c["time"]["duration_s"]) + dt_box, dt_box)
    efac, gfac, kdil = (np.interp(st, times, met[k]) for k in (5, 6, 7))
    tl = " ".join(f"{t:.1f}" for t in st)
    (d / "aero_emit.dat").write_text(f"time {tl}\nrate {' '.join(f'{x:.6e}' for x in efac)}\n"
                                     f"dist {' '.join(['aero_emit_dist.dat'] * len(st))}\n")
    (d / "aero_back.dat").write_text(f"time {tl}\nrate {' '.join(f'{x:.6e}' for x in kdil)}\n"
                                     f"dist {' '.join(['aero_back_dist.dat'] * len(st))}\n")
    ge = c["emissions"].get("gas_ugm3_s") or {}
    lines = [f"time {tl}", f"rate {' '.join(f'{x:.6e}' for x in gfac)}"]
    for s in gases:     # mol m-2 s-1 with a 1 m column = mol m-3 s-1
        lines.append(f"{GAS_MAP[s]} " + " ".join(f"{float(ge.get(s, 0.0)) / (GAS_MW[s] * 1e6):.6e}" for _ in st))
    (d / "gas_emit.dat").write_text("\n".join(lines) + "\n")
    lines = [f"time {tl}", f"rate {' '.join(f'{x:.6e}' for x in kdil)}"]
    for s in gases:
        lines.append(f"{GAS_MAP[s]} " + " ".join(f"{_ppb(float(gb.get(s, 0.0)), GAS_MW[s], T0, P0):.6e}" for _ in st))
    (d / "gas_back.dat").write_text("\n".join(lines) + "\n")
    tl = " ".join(f"{t:.1f}" for t in times)
    (d / "temp.dat").write_text(f"time {tl}\ntemp {' '.join(f'{x:.4f}' for x in met[1])}\n")
    (d / "pres.dat").write_text(f"time {tl}\npressure {' '.join(f'{x:.2f}' for x in met[3])}\n")
    (d / "height.dat").write_text("time 0\nheight 1\n")
    _series_file(d / "rh.dat", times, met[2])
    _series_file(d / "jno2.dat", times, met[4] if c["processes"].get("soa_photolysis", True) else 0 * met[4])

    p = c["processes"]
    if p["nucleation"]["enabled"]:
        print("warning: nucleation is not represented in the PartMC reference run")
    if (c.get("chemistry") or {}).get("n2o5_uptake", {}).get("enabled"):
        print("warning: online N2O5 uptake is not represented in the PartMC reference run")
    t = c["time"]
    out_every = float(t["output_every_s"])
    spec = f"""run_type particle
output_prefix out/pmc
n_repeat 1
n_part {int(n_part)}
restart no
do_select_weighting yes
weight_type power
weighting_exponent {weighting_exponent}

t_max {float(t['duration_s'])}
del_t {float(del_t)}
t_output {out_every}
t_progress {max(out_every * 4, 3600.0)}

do_camp_chem no
do_tchem no

gas_data gas_data.dat
gas_init gas_init.dat

aerosol_data aero_data.dat
do_fractal no
aerosol_init aero_init_dist.dat

temp_profile temp.dat
pressure_profile pres.dat
height_profile height.dat
gas_emissions gas_emit.dat
gas_background gas_back.dat
aero_emissions aero_emit.dat
aero_background aero_back.dat
loss_function none

rel_humidity {float(met[2][0])}
latitude 0
longitude 0
altitude 0
start_time 0
start_day 200

do_coagulation {'yes' if p['coagulation'] else 'no'}
{'coag_kernel brown' if p['coagulation'] else ''}
do_condensation no
do_mosaic no
do_nucleation no
do_immersion_freezing no

rand_init {int(seed)}
allow_doubling yes
allow_halving yes
record_removals no
do_parallel no
"""
    (d / "pmc.spec").write_text("\n".join(x for x in spec.splitlines() if x.strip() != "") + "\n")
    soa = any(s.startswith("PSOA") for s in table.aerosol)
    (d / "thermo.nml").write_text(
        "&scram_thermo\n  enabled = %s\n  courant = 0.05\n  alfhp = 0.1\n  soa = %s\n  crustal = .true.\n"
        "  dust_species = 'OIN'\n  jno2_file = 'jno2.dat'\n  rh_file = 'rh.dat'\n/\n"
        % (".true." if p["condensation"] else ".false.", ".true." if soa else ".false."))
    casemod.save(c, d / "case.resolved.json")
    return c


def run_partmc(user_case: dict, run_dir, n_part: int = 10000, del_t: float = 60.0, seed: int = 1,
               timeout_s: float = 6 * 3600) -> dict:
    d = Path(run_dir).resolve()
    c = write_partmc(user_case, d, n_part=n_part, del_t=del_t, seed=seed)
    (d / "out").mkdir(exist_ok=True)
    exe = PARTMC_BIN["v23" if c["kernel"]["isorropia"] == "v23" else "v1"]
    env = dict(os.environ, PMC_SCRAM_THERMO_NML=str(d / "thermo.nml"))
    t0 = time.time()
    with open(d / "pmc.log", "w") as log:
        proc = subprocess.run([exe, "pmc.spec"], cwd=d, env=env, stdout=log, stderr=subprocess.STDOUT,
                              timeout=timeout_s)
    wall = time.time() - t0
    if proc.returncode != 0:
        raise RuntimeError("PartMC failed: " + (d / "pmc.log").read_text()[-1500:])
    info = project(d)
    info["summary"]["wall_s"] = wall
    info["summary"]["n_part"] = n_part
    (d / "summary.json").write_text(json.dumps(info["summary"], indent=1))
    truth(d)
    return info


# ------------------------------------------------------------------ reading PartMC output
def read_particles(run_dir):
    """Yield (time, env dict, species names, mass (n_spec, n_part) kg, num_conc (n_part) m-3,
    gas names, gas ppb) for every PartMC output file."""
    from netCDF4 import Dataset
    for f in sorted(glob.glob(str(Path(run_dir) / "out" / "pmc_0001_*.nc"))):
        nc = Dataset(f)
        spec = nc.variables["aero_species"].names.split(",")
        gas_names = nc.variables["gas_species"].names.split(",")
        env = {k: float(nc.variables[k][...]) for k in ("temperature", "relative_humidity", "pressure")
               if k in nc.variables}
        yield (float(nc.variables["time"][0]), env, spec, np.array(nc.variables["aero_particle_mass"][:]),
               np.array(nc.variables["aero_num_conc"][:]), gas_names, np.array(nc.variables["gas_mixing_ratio"][:]))
        nc.close()


def project(run_dir, out_dir=None, mixing: dict | None = None) -> dict:
    """PartMC particles -> SCRAM bin x class arrays, stored like a SCRAM run (in run_dir, or
    in out_dir with the case's mixing replaced by `mixing`, e.g. {"preset": "NC11"})."""
    d = Path(run_dir)
    c = casemod.load(d / "case.resolved.json")
    if mixing is not None:
        c = casemod.set_path(c, "mixing", mixing)
    od = Path(out_dir) if out_dir else d
    od.mkdir(parents=True, exist_ok=True)
    g = build_grid(c)
    table, bounds, classes = g["table"], np.array(g["bounds_um"]), g["classes"]
    nsize, nfrac, naer = len(bounds) - 1, len(classes), len(table.aerosol)
    rh_t = np.loadtxt(d / "rh.dat")
    times, met, NUM, MASS, GAS, OUT = [], [], [], [], [], []
    inv = {_pmc_name(s): i for i, s in enumerate(table.aerosol)}
    gidx = {GAS_MAP[s]: table.gases.index(s) for s in GAS_MAP if s in table.gases}
    dens = np.array([table.density_gcm3(s) for s in table.aerosol]) * 1e3      # kg m-3
    iw = table.aerosol.index("PH2O")
    group_of = [table.group_index(s) for s in table.aerosol]                   # 0 = water
    for t, env, spec, mass, nconc, gnames, gppb in read_particles(d):
        m = np.zeros((naer, mass.shape[1]))
        for k, s in enumerate(spec):
            if s in inv:
                m[inv[s]] = mass[k]
        dry = m.sum(axis=0) - m[iw]
        vdry = (np.delete(m, iw, axis=0) / np.delete(dens, iw)[:, None]).sum(axis=0)
        ddry_um = (6 * vdry / math.pi) ** (1 / 3) * 1e6
        b = np.clip(np.searchsorted(bounds, ddry_um, side="right") - 1, 0, nsize - 1)
        out_of_range = float(nconc[(ddry_um < bounds[0]) | (ddry_um > bounds[-1])].sum())
        num = np.zeros((nsize, nfrac))
        ms = np.zeros((nsize, nfrac, naer))
        for p in range(m.shape[1]):
            if dry[p] <= 0:
                continue
            gf = np.zeros(len(table.groups))
            for i in range(naer):
                if group_of[i] > 0:
                    gf[group_of[i] - 1] += m[i, p]
            cls = find_class(list(gf / dry[p]), classes)
            num[b[p], cls] += nconc[p]
            ms[b[p], cls] += m[:, p] * nconc[p] * 1e9                           # ug m-3
        T, P = env.get("temperature", 298.0), env.get("pressure", 1e5)
        rh = float(np.interp(t, rh_t[:, 0], rh_t[:, 1]))
        gas = np.zeros(len(table.gases))
        for gn, val in zip(gnames, gppb):
            if gn in gidx:
                s = table.gases[gidx[gn]]
                gas[gidx[gn]] = val * 1e-9 * P / (R * T) * GAS_MW[s] * 1e6
        times.append(t); met.append([T, rh, P]); NUM.append(num); MASS.append(ms); GAS.append(gas)
        OUT.append(out_of_range)
    res = {"time_s": np.array(times), "met": np.array(met), "number_m3": np.array(NUM),
           "mass_ugm3": np.array(MASS), "gas_ugm3": np.array(GAS)}
    d = od
    np.savez_compressed(d / "results.npz", **res)
    casemod.save(c, d / "case.resolved.json")
    from .grid import class_label
    gridinfo = {"bounds_um": list(bounds), "nsize": nsize, "nfrac": nfrac,
                "classes": [[list(x) for x in cls] for cls in classes],
                "class_labels": [class_label(cls, table.groups) for cls in classes],
                "species": table.to_dict(), "kernel": "partmc-" + c["kernel"]["isorropia"],
                "notes": [f"PartMC particles projected on the case grid; number outside the size range at the end: "
                          f"{OUT[-1] / 1e6:.3g} cm-3" if OUT else ""]}
    (d / "grid.json").write_text(json.dumps(gridinfo, indent=1))
    summary = {"name": c["case"]["name"] + " (PartMC)", "returncode": 0, "kernel": "PartMC + pmc_scram_thermo",
               "cpu_s": None, "n_output": len(times), "notes": gridinfo["notes"]}
    (d / "summary.json").write_text(json.dumps(summary, indent=1))
    return {"case": c, "grid": gridinfo, "results": res, "summary": summary}


def truth(run_dir, supersats=(0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0)) -> dict:
    """Particle-resolved metrics (no binning): per-particle kappa-Koehler CCN spectrum
    (same Koehler form and kappa table as CESM-SCRAM), mixing-state index chi over the
    dry species (Riemer and West 2013), BC-containing particle number and coating ratio."""
    d = Path(run_dir)
    c = casemod.load(d / "case.resolved.json")
    table = build_grid(c)["table"]
    dens = {_pmc_name(s): table.density_gcm3(s) * 1e3 for s in table.aerosol}
    kap = {_pmc_name(s): KAPPA.get(s, 0.1) for s in table.aerosol}
    out = {"time_s": [], "ccn_cm3": {str(s): [] for s in supersats}, "chi": [], "D_alpha": [], "D_gamma": [],
           "n_bc_cm3": [], "bc_coating_ratio": [], "n_total_cm3": []}
    for t, env, spec, mass, nconc, *_ in read_particles(d):
        T = env.get("temperature", 298.0)
        iw = spec.index("H2O")
        dry_idx = [k for k in range(len(spec)) if k != iw]
        md = mass[dry_idx]                               # kg
        mtot = md.sum(axis=0)
        ok = mtot > 0
        md, nc_, mtot = md[:, ok], nconc[ok], mtot[ok]
        names = [spec[k] for k in dry_idx]
        vol = np.array([md[i] / dens[n] for i, n in enumerate(names)])
        vtot = vol.sum(axis=0)
        keff = (vol * np.array([kap[n] for n in names])[:, None]).sum(axis=0) / vtot
        r = (3 * vtot / (4 * math.pi)) ** (1 / 3)
        aten = 2 * 18.016 * 0.076 / (8.31446261815324e3 * T * 1000.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            sc = np.where(keff > 0, 2 * aten * np.sqrt(aten / (27 * keff * r ** 3)), np.inf)
        for s in supersats:
            out["ccn_cm3"][str(s)].append(float(nc_[sc <= s / 100.0].sum() / 1e6))
        # mixing-state index over dry species, number-weighted by num_conc x particle mass
        w = nc_ * mtot
        p_ia = md / mtot                                 # species mass fractions in each particle
        with np.errstate(divide="ignore", invalid="ignore"):
            H_i = -np.nansum(np.where(p_ia > 0, p_ia * np.log(p_ia), 0.0), axis=0)
        p_i = w / w.sum()
        D_alpha = float(np.exp((p_i * H_i).sum()))
        p_a = (md * nc_).sum(axis=1) / (md * nc_).sum()
        D_gamma = float(np.exp(-np.sum(np.where(p_a > 0, p_a * np.log(np.where(p_a > 0, p_a, 1.0)), 0.0))))
        chi = (D_alpha - 1) / (D_gamma - 1) if D_gamma > 1 else 1.0
        bc = names.index("BC") if "BC" in names else None
        if bc is not None:
            has = md[bc] > 0
            out["n_bc_cm3"].append(float(nc_[has].sum() / 1e6))
            bcm = (md[bc] * nc_).sum()
            out["bc_coating_ratio"].append(float(((mtot - md[bc]) * nc_)[has].sum() / bcm) if bcm > 0 else None)
        out["time_s"].append(t)
        out["chi"].append(chi); out["D_alpha"].append(D_alpha); out["D_gamma"].append(D_gamma)
        out["n_total_cm3"].append(float(nc_.sum() / 1e6))
    (d / "truth.json").write_text(json.dumps(out))
    return out


def chi_sections(mass_ugm3: np.ndarray, number_m3: np.ndarray, species: list[str]) -> float:
    """Mixing-state index chi of a sectional state, each (bin, class) section treated as
    one particle type with its mean dry composition. With identical inputs this is the
    mixing state a representation can hold (NC1 keeps only size-dependent diversity)."""
    iw = species.index("PH2O")
    m = np.delete(mass_ugm3.reshape(-1, mass_ugm3.shape[-1]), iw, axis=1)   # (sections, species)
    mt = m.sum(axis=1)
    ok = mt > 0
    m, mt = m[ok], mt[ok]
    p = m / mt[:, None]
    with np.errstate(divide="ignore", invalid="ignore"):
        H = -np.where(p > 0, p * np.log(p), 0.0).sum(axis=1)
    w = mt / mt.sum()
    D_alpha = float(np.exp((w * H).sum()))
    pa = m.sum(axis=0) / m.sum()
    D_gamma = float(np.exp(-np.where(pa > 0, pa * np.log(np.where(pa > 0, pa, 1.0)), 0.0).sum()))
    return (D_alpha - 1) / (D_gamma - 1) if D_gamma > 1 else 1.0
