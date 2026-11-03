"""Compile a case.json into the inputs of the Fortran driver (kernel/driver/scram_box.F90).

Writes into a run directory:
  case.resolved.json   the case with all defaults filled in
  grid.json            sections, classes, species and gases (used to read the results)
  scram_species.dat    species table for the kernel (and scram_species_noiso.dat)
  run.nml              scalars and switches
  inputs.txt           arrays as keyword blocks
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from . import case as casemod
from . import grid
from . import species as sp

UGUM3_PER_GCM3 = 1.0e-6          # SCRAM density unit: ug/um3
REDISTRIBUTION = {n: i for i, n in enumerate(casemod.ENUMS["processes.redistribution"])}
SOLVER = {"euler": 0, "etr": 1, "ros2": 2}


def species_table(c: dict) -> sp.SpeciesTable:
    profile = "noiso" if c["kernel"]["isorropia"] == "none" else "iso"
    t = sp.with_custom_species(sp.library(profile), c["species"].get("custom", []))
    g = c["species"].get("groups")
    if g:
        t = sp.with_groups(t, g["groups"], g["group_of"])
    return t


def build_grid(c: dict) -> dict:
    table = species_table(c)
    bounds = grid.size_bounds_um(c["size_bins"])
    classes = grid.composition_classes(len(table.groups), c["mixing"], list(table.groups))
    return {"table": table, "bounds_um": bounds, "classes": classes}


def _dry_density(fractions: dict, table: sp.SpeciesTable) -> float:
    inv = sum(f / table.density_gcm3(s) for s, f in fractions.items())
    return 1.0 / inv


def _normalised(fractions: dict, table: sp.SpeciesTable, where: str) -> dict:
    fr = {k: float(v) for k, v in (fractions or {}).items() if float(v) > 0}
    for s in fr:
        if s not in table.aerosol or s == "PH2O":
            raise ValueError(f"{where}: unknown or non-dry species '{s}'")
    tot = sum(fr.values())
    if tot <= 0:
        raise ValueError(f"{where}: mass_fraction is empty")
    return {k: v / tot for k, v in fr.items()}


def _group_fractions(fr: dict, table: sp.SpeciesTable) -> list[float]:
    g = [0.0] * len(table.groups)
    for s, f in fr.items():
        g[table.group_index(s) - 1] += f
    return g


def particles_to_sections(spec: dict, g: dict, rate: bool, where: str):
    """Map log-normal modes (or a class x bin table) to section number (#/m3[/s]) and
    mass (ug/m3[/s]) arrays shaped (nbin, naer); returns (number, mass, notes)."""
    table, bounds, classes = g["table"], g["bounds_um"], g["classes"]
    nsize, nfrac, naer = len(bounds) - 1, len(classes), len(table.aerosol)
    number = np.zeros(nsize * nfrac)
    mass = np.zeros((nsize * nfrac, naer))
    notes = []
    key_n = "number_cm3_s" if rate else "number_cm3"
    mode = spec.get("mode", "lognormal")
    if mode == "lognormal":
        for k, m in enumerate(spec.get("modes", [])):
            fr = _normalised(m.get("mass_fraction"), table, f"{where}.modes[{k}]")
            rho = _dry_density(fr, table)                      # g/cm3
            cls = grid.find_class(_group_fractions(fr, table), classes)
            nb, vb, out_n, out_v = grid.lognormal_to_bins(float(m[key_n]) * 1e6, float(m["dg_um"]),
                                                          float(m["sigma_g"]), bounds)
            if out_v > 0.01:
                notes.append(f"{where}.modes[{k}] '{m.get('name', k)}': {100*out_v:.1f}% of the "
                             f"volume lies outside the size range and is dropped")
            for b in range(nsize):
                j = b * nfrac + cls
                number[j] += nb[b]
                mdry = vb[b] * rho * UGUM3_PER_GCM3                 # ug/m3
                for s, f in fr.items():
                    mass[j, table.aerosol.index(s)] += mdry * f
    elif mode == "bin_mass":
        # mass per species and size bin (ug/m3[/s]); number per bin optional (else from
        # the dry volume at the geometric mid diameter). "mixing": "internal" puts each
        # bin into the class of its total composition; "by_group" makes one population
        # per species group, "by_group_set" one per listed species set ("sets"), each in
        # the class of its own composition with the number split by volume (the legacy
        # SCRAM external start).
        key_m = "mass_ugm3_s" if rate else "mass_ugm3"
        masses = {s: np.asarray(v, dtype=float) for s, v in spec[key_m].items()}
        for s, v in masses.items():
            if s not in table.aerosol or s == "PH2O":
                raise ValueError(f"{where}.{key_m}: unknown or non-dry species '{s}'")
            if v.shape != (nsize,):
                raise ValueError(f"{where}.{key_m}.{s} needs {nsize} values")
        nin = spec.get(key_n)
        nin = None if nin is None else np.asarray(nin, dtype=float) * 1e6
        how = spec.get("mixing", "internal")
        if how == "internal":
            sets = [list(masses)]
        elif how == "by_group":
            sets = [[s for s in masses if table.group_of.get(s) == gname] for gname in table.groups]
        elif how == "by_group_set":
            sets = [[s for s in st if s in masses] for st in spec["sets"]]
            missing = set(masses) - {s for st in sets for s in st}
            if missing:
                raise ValueError(f"{where}.sets does not cover {sorted(missing)}")
        else:
            raise ValueError(f"{where}.mixing must be 'internal', 'by_group' or 'by_group_set'")
        for b in range(nsize):
            dmid = math.sqrt(bounds[b] * bounds[b + 1])
            vol = {s: v[b] / (table.density_gcm3(s) * UGUM3_PER_GCM3) for s, v in masses.items()}
            vtot = sum(vol.values())
            if vtot <= 0:
                continue
            ntot = nin[b] if nin is not None else vtot / (math.pi / 6 * dmid ** 3)
            parts = []
            for members in sets:
                vs = sum(vol[s] for s in members)
                ms = sum(masses[s][b] for s in members)
                if vs > 0:
                    fr = {s: masses[s][b] / ms for s in members}
                    parts.append((grid.find_class(_group_fractions(fr, table), classes), members, vs / vtot))
            for cls, members, share in parts:
                j = b * nfrac + cls
                number[j] += ntot * share
                for s in members:
                    mass[j, table.aerosol.index(s)] += masses[s][b]
    elif mode == "table":
        n = np.asarray(spec[key_n], dtype=float)               # (nclass, nsize) #/cm3
        if n.shape != (nfrac, nsize):
            raise ValueError(f"{where}.{key_n} must be {nfrac} classes x {nsize} bins")
        comp = spec["mass_fraction"]                           # list per class of dicts
        if len(comp) != nfrac:
            raise ValueError(f"{where}.mass_fraction needs one entry per class")
        for c_i in range(nfrac):
            fr = _normalised(comp[c_i], table, f"{where}.mass_fraction[{c_i}]")
            rho = _dry_density(fr, table)
            for b in range(nsize):
                j = b * nfrac + c_i
                dmid = math.sqrt(bounds[b] * bounds[b + 1])     # um
                number[j] += n[c_i, b] * 1e6
                mdry = n[c_i, b] * 1e6 * math.pi / 6 * dmid ** 3 * rho * UGUM3_PER_GCM3
                for s, f in fr.items():
                    mass[j, table.aerosol.index(s)] += mdry * f
    else:
        raise ValueError(f"{where}.mode must be 'lognormal', 'bin_mass' or 'table'")
    return number, mass, notes


def gas_vector(values: dict, table: sp.SpeciesTable, where: str) -> np.ndarray:
    v = np.zeros(len(table.gases))
    for k, x in (values or {}).items():
        if k not in table.gases:
            raise ValueError(f"{where}: unknown gas '{k}'")
        v[table.gases.index(k)] = float(x)
    return v


def met_table(c: dict) -> np.ndarray:
    """Rows: t, T, RH, P, J_NO2, particle emission factor, gas source factor,
    dilution rate (1/s), on the union of all time axes."""
    e, em = c["environment"], c["emissions"]
    prof, gprof = em.get("time_profile"), em.get("gas_time_profile")
    ts = e.get("time_series") or {}
    axes = [0.0] + [float(x) for x in ts.get("time_s", [])]
    for p in (prof, gprof):
        if isinstance(p, dict):
            axes += [float(x) for x in p["time_s"]]
    times = sorted(set(axes))

    def series(name, const):
        if ts.get(name) is not None and ts.get("time_s"):
            return lambda t: float(np.interp(t, ts["time_s"], ts[name]))
        return lambda t: float(const)

    def profile(p):
        if isinstance(p, dict):
            return lambda t: float(np.interp(t, p["time_s"], p["factor"]))
        return lambda t: 1.0

    f_t = series("temperature_K", e["temperature_K"])
    f_rh = series("relative_humidity", e["relative_humidity"])
    f_p = series("pressure_Pa", e["pressure_Pa"])
    f_j = series("jno2_s1", e.get("jno2_s1", 0.0))
    f_k = series("dilution_s1", (e.get("dilution") or {}).get("rate_s1", 0.0))
    het = (c.get("chemistry") or {}).get("n2o5_uptake") or {}
    if not het.get("enabled"):
        f_n = lambda t: 0.0
    elif het.get("mode", "production") == "production":
        f_n = series("n2o5_production_ugm3_s", het.get("n2o5_production_ugm3_s", 0.0))
    else:
        f_n = series("n2o5_ugm3", het.get("n2o5_ugm3", 0.0))
    f_e = profile(prof)
    if gprof == "jno2":                       # photochemical production proxy
        jmax = max([f_j(t) for t in times] + [1e-30])
        f_g = lambda t: f_j(t) / jmax
    elif gprof is None:
        f_g = f_e
    else:
        f_g = profile(gprof)
    return np.array([[t, f_t(t), min(max(f_rh(t), 0.0), 1.0), f_p(t), max(f_j(t), 0.0),
                      f_e(t), f_g(t), max(f_k(t), 0.0), max(f_n(t), 0.0)] for t in times]).T


CHM_SPECIES = {"PSO4", "PNH4", "PNO3", "PBC", "PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P"}
MAXSP = 64                                   # size of the driver's per-species het arrays


def het_namelist(c: dict, table: sp.SpeciesTable) -> dict:
    """&box_het for the driver's online N2O5 uptake."""
    h = (c.get("chemistry") or {}).get("n2o5_uptake") or {}
    if len(table.aerosol) > MAXSP:
        raise ValueError(f"at most {MAXSP} aerosol species")
    sg = h.get("species_gamma")
    if sg is None:                                   # the N2O5 column of the post-processing table
        sg = (((c.get("postprocess") or {}).get("hetchem") or {}).get("species_gamma") or {}).get("N2O5")
    if sg is None:                                   # default reproduces the CESM selection
        sg = {s: (float(h.get("gamma_cesm", 0.1)) if s in CHM_SPECIES else 0.0) for s in table.aerosol}
    return {"het_n2o5": int(bool(h.get("enabled"))),
            "het_rule": 0 if h.get("rule", "cesm") == "cesm" else 1,
            "clno2_mode": 1 if h.get("clno2_yield", "bertram_thornton") == "bertram_thornton" else 0,
            "gamma_cesm": float(h.get("gamma_cesm", 0.1)),
            "gamma_n2o5": [float(sg.get(s, 0.0)) for s in table.aerosol],
            "is_chm": [s in CHM_SPECIES for s in table.aerosol],
            "n2o5_mode": 1 if h.get("mode", "production") == "production" else 0,
            "k_other_s1": float(h.get("k_other_s1", 0.0))}


def _block(name, arr) -> str:
    a = np.asarray(arr).ravel(order="F")
    if a.dtype.kind in "iu":
        body = " ".join(str(int(x)) for x in a)
    else:
        body = "\n".join(" ".join(f"{x:.10e}" for x in a[i:i + 5]) for i in range(0, len(a), 5))
    return f"{name} {a.size}\n{body}\n"


def compile_case(user_case: dict, run_dir: str | Path) -> dict:
    c = casemod.resolve(user_case)
    run = Path(run_dir)
    run.mkdir(parents=True, exist_ok=True)
    g = build_grid(c)
    table, bounds, classes = g["table"], g["bounds_um"], g["classes"]
    nsize, nfrac, naer, ngas = len(bounds) - 1, len(classes), len(table.aerosol), len(table.gases)

    n0, m0, notes0 = particles_to_sections(c["initial"]["particles"], g, False, "initial.particles")
    ne, me, notes1 = particles_to_sections(c["emissions"]["particles"], g, True, "emissions.particles")
    gas0 = gas_vector(c["initial"].get("gas_ugm3"), table, "initial.gas_ugm3")
    gsrc = gas_vector(c["emissions"].get("gas_ugm3_s"), table, "emissions.gas_ugm3_s")
    met = met_table(c)
    bg = (c["environment"].get("dilution") or {}).get("background", "initial")
    if bg in (None, "initial"):
        nb, mb, gb = n0, m0, gas0
    else:
        nb, mb, notes2 = particles_to_sections(bg.get("particles") or {"modes": []}, g, False,
                                               "environment.dilution.background.particles")
        notes0 = notes0 + notes2
        gb = gas_vector(bg.get("gas_ugm3"), table, "environment.dilution.background.gas_ugm3")

    # species tables: the wrapper appends _noiso to the file name when ISORROPIA is off
    sp.write_species_file(table, run / "scram_species.dat")
    if c["kernel"]["isorropia"] == "none":
        sp.write_species_file(table, run / "scram_species_noiso.dat")

    ids = {s: (table.aerosol.index(s) + 1 if s in table.aerosol else -1)
           for s in ["PMD", "PBC", "PNA", "PSO4", "PNH4", "PNO3", "PHCL", "PPOAlP", "PH2O",
                     "PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P"]}
    p, a, t = c["processes"], c["advanced"], c["time"]
    nml = {
        "box_dims": {"nsize": nsize, "nfrac": nfrac, "naer": naer, "ngas": ngas,
                     "ngroup": len(table.groups)},
        "box_time": {"duration_s": t["duration_s"], "dt_s": t["dt_s"], "substep_s": t["substep_s"],
                     "output_every_s": t["output_every_s"]},
        "box_proc": {
            "with_coagulation": int(bool(p["coagulation"])),
            "with_condensation": int(bool(p["condensation"])),
            "with_nucleation": int(bool(p["nucleation"]["enabled"])),
            "nucleation_model": 1 if p["nucleation"]["model"] == "ternary" else 0,
            "dynamic_solver": SOLVER[p["dynamic_solver"]],
            "thermodynamic_model": 0 if p["thermodynamics"] == "dynamic" else 1,
            "use_isorropia": 0 if c["kernel"]["isorropia"] == "none" else 1,
            "sulfate_computation": 1 if p["sulfate"] == "dynamic" else 0,
            "redistribution_method": int(p["redistribution_method"]),
            "rdb_core_conserv": REDISTRIBUTION[p["redistribution"]],
            "wet_diam_estimation": 0 if p["wet_diameter"] == "isorropia" else 1,
            "density_mode": 0 if p["density"]["mode"] == "variable" else 1,
            "fixed_density_kgm3": p["density"]["fixed_kgm3"],
            "bc_core_floor_um": a["bc_core_floor_um"], "salt_core_floor_um": a["salt_core_floor_um"],
            "dust_core_floor_um": a["dust_core_floor_um"],
            "salt_rdb_bidirectional": int(bool(a["salt_rdb_bidirectional"])),
            "salt_fix_size_rdb": int(bool(a["salt_fix_size_rdb"])),
            "dust_fix_size_rdb": int(bool(a["dust_fix_size_rdb"])),
            "nucl_rate_cap_cm3_s": a["nucl_rate_cap_cm3_s"],
            "nucl_h2so4_gas_fraction_max": a["nucl_h2so4_gas_fraction_max"],
            "nucl_nh3_gas_fraction_max": a["nucl_nh3_gas_fraction_max"],
            "soa_photolysis": int(bool(p.get("soa_photolysis", True))),
            # written only when used, so kernels built before the option still read the namelist
            # split: assumed composition spread alpha (fraction of the class width), in per mille
            **({"class_split": max(1, round(1000 * min(1.0, max(0.001, float(p.get("class_split_alpha", 1.0))))))}
               if p.get("class_redistribution", "whole") == "split" else {}),
        },
        "box_ids": {"sid_pmd": ids["PMD"], "sid_pbc": ids["PBC"], "sid_pna": ids["PNA"],
                    "sid_pso4": ids["PSO4"], "sid_pnh4": ids["PNH4"], "sid_pno3": ids["PNO3"],
                    "sid_phcl": ids["PHCL"], "sid_poa": ids["PPOAlP"], "sid_ph2o": ids["PH2O"],
                    "sid_soa": [ids[f"PSOA{i}P"] for i in range(1, 6)],
                    "gid_sulf": table.gases.index("SULF") + 1 if "SULF" in table.gases else -1,
                    "gid_hno3": table.gases.index("HNO3") + 1 if "HNO3" in table.gases else -1},
        "box_het": het_namelist(c, table),
        "box_files": {"inputs_file": "inputs.txt", "output_file": "output.txt"},
    }
    lines = []
    fval = lambda x: (".true." if x else ".false.") if isinstance(x, bool) else (f"{x:.10g}" if isinstance(x, float) else str(x))
    for group, vals in nml.items():
        lines.append(f"&{group}")
        for k, v in vals.items():
            if isinstance(v, list):
                v = ", ".join(fval(x) for x in v)
            elif isinstance(v, bool):
                v = fval(v)
            elif isinstance(v, str):
                v = f"'{v}'"
            elif isinstance(v, float):
                v = f"{v:.10g}"
            lines.append(f"  {k} = {v}")
        lines.append("/")
    (run / "run.nml").write_text("\n".join(lines) + "\n")

    density = [table.density_gcm3(s) * UGUM3_PER_GCM3 for s in table.aerosol]
    group_rel = np.array([table.group_index(s) for s in table.aerosol], dtype=int)
    with open(run / "inputs.txt", "w") as f:
        f.write(_block("bin_bound_m", np.array(bounds) * 1e-6))
        f.write(_block("comp_bounds", grid.comp_bounds_array(classes)))
        f.write(_block("density_ugum3", density))
        f.write(_block("group_rel", group_rel))
        f.write(_block("number_m3", n0))
        f.write(_block("mass_ugm3", m0))
        f.write(_block("gas_ugm3", gas0))
        f.write(_block("emis_number_m3s", ne))
        f.write(_block("emis_mass_ugm3s", me))
        f.write(_block("gas_source_ugm3s", gsrc))
        f.write(_block("bg_number_m3", nb))
        f.write(_block("bg_mass_ugm3", mb))
        f.write(_block("bg_gas_ugm3", gb))
        f.write(f"met {met.shape[0]} {met.shape[1]}\n")
        for col in met.T:
            f.write(" ".join(f"{x:.10e}" for x in col) + "\n")

    gridinfo = {
        "bounds_um": bounds, "nsize": nsize, "nfrac": nfrac,
        "classes": [[list(b) for b in cls] for cls in classes],
        "class_labels": [grid.class_label(cls, table.groups) for cls in classes],
        "species": table.to_dict(), "kernel": c["kernel"]["isorropia"],
        "notes": notes0 + notes1,
    }
    (run / "grid.json").write_text(json.dumps(gridinfo, indent=1, ensure_ascii=False))
    casemod.save(c, run / "case.resolved.json")
    return {"case": c, "grid": gridinfo, "run_dir": str(run)}
