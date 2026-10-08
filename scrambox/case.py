"""case.json: defaults, merging, validation and variant paths."""
from __future__ import annotations

import copy
import json
from pathlib import Path

SCHEMA = "scram-box/0.1"

# Defaults follow the CESM-SCRAM production profile (NC6 cfg) except where a box
# study needs otherwise (ISORROPIA II 2.3 with the dynamic scheme).
DEFAULTS = {
    "schema": SCHEMA,
    "case": {"name": "untitled", "description": "", "scenario_id": None, "category": "custom",
             "reference": "", "derived_from": None},
    "kernel": {"isorropia": "v23", "crustal_speciation": "cmaq_aero6_7"},
    "environment": {"temperature_K": 298.0, "pressure_Pa": 101325.0, "relative_humidity": 0.70,
                    "jno2_s1": 0.0, "time_series": None,
                    "dilution": {"rate_s1": 0.0, "background": "initial"}},
    "time": {"duration_s": 43200.0, "dt_s": 600.0, "substep_s": None, "output_every_s": 600.0},
    "species": {"custom": []},
    "mixing": {"preset": "NC6", "external": True,
               "primary_fraction_bounds": [0.0, 0.5, 1.0],
               "secondary_fraction_bounds": [0.0, 0.5, 1.0]},
    "size_bins": {"mode": "explicit",
                  "bounds_um": [0.003, 0.01, 0.0398, 0.1585, 0.631, 2.5, 5.0, 10.0]},
    "initial": {"particles": {"mode": "lognormal", "modes": []}, "gas_ugm3": {}},
    "emissions": {"particles": {"mode": "lognormal", "modes": []}, "gas_ugm3_s": {},
                  "time_profile": None, "gas_time_profile": None},
    "processes": {"coagulation": True, "condensation": True,
                  "nucleation": {"enabled": False, "model": "binary"},
                  "sulfate": "dynamic", "thermodynamics": "dynamic", "dynamic_solver": "etr",
                  "redistribution_method": 3, "redistribution": "moving_center_dualpivot",
                  "wet_diameter": "isorropia", "density": {"mode": "variable", "fixed_kgm3": 1400.0},
                  "soa_photolysis": True},
    "advanced": {"bc_core_floor_um": 0.05, "salt_core_floor_um": 0.05, "dust_core_floor_um": 0.90,
                 "salt_rdb_bidirectional": False, "salt_fix_size_rdb": False,
                 "dust_fix_size_rdb": False, "nucl_rate_cap_cm3_s": 1.0e3,
                 "nucl_h2so4_gas_fraction_max": 0.25, "nucl_nh3_gas_fraction_max": 0.25},
    # online heterogeneous N2O5 uptake; mode "production": P(N2O5) with steady state
    # [N2O5] = P / (k_het + k_other); mode "concentration": prescribed [N2O5].
    # environment.time_series may carry "n2o5_production_ugm3_s" or "n2o5_ugm3".
    # Products (2 - phi) HNO3 and phi ClNO2 (particle chloride).
    "chemistry": {"n2o5_uptake": {"enabled": False, "mode": "production",
                                  "n2o5_production_ugm3_s": 1.2e-3, "k_other_s1": 2.0e-3,
                                  "n2o5_ugm3": 0.0, "rule": "cesm",
                                  "gamma_cesm": 0.1, "species_gamma": None,
                                  "clno2_yield": "bertram_thornton"}},
    "postprocess": {"optics": {"enabled": False}, "ccn": {"enabled": False},
                    "hetchem": {"enabled": False}},
    "output": {"format": "npz"},
}

ENUMS = {
    "kernel.isorropia": ["v23", "v1", "none"],
    "processes.sulfate": ["dynamic", "equilibrium"],
    "processes.thermodynamics": ["dynamic", "equilibrium"],
    "processes.dynamic_solver": ["euler", "etr", "ros2"],
    "processes.redistribution": ["legacy", "core_conserv", "core_nogrow", "core_smallgrow",
                                 "core_fluxform", "core_ppm", "core_ppm_dualmoment",
                                 "core_ppm_cutnum", "moving_center_dualpivot"],
    "processes.wet_diameter": ["isorropia", "gerber"],
    "processes.nucleation.model": ["binary", "ternary"],
    "mixing.preset": ["NC1", "NC6", "NC11", "custom"],
    "chemistry.n2o5_uptake.rule": ["cesm", "species"],
    "chemistry.n2o5_uptake.mode": ["production", "concentration"],
    "chemistry.n2o5_uptake.clno2_yield": ["none", "bertram_thornton"],
}


def deep_merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def get_path(d: dict, path: str):
    cur = d
    for p in path.split("."):
        cur = cur[p]
    return cur


def set_path(d: dict, path: str, value) -> dict:
    out = copy.deepcopy(d)
    cur = out
    parts = path.split(".")
    for p in parts[:-1]:
        cur = cur.setdefault(p, {})
    cur[parts[-1]] = copy.deepcopy(value)
    return out


def substep_default(case: dict) -> float:
    """Kernel sub-step: 10 s for the dynamic (non-equilibrium) ISORROPIA scheme,
    600 s (CESM scram_target_dt) for no-ISORROPIA and bulk-equilibrium runs."""
    if case["kernel"]["isorropia"] != "none" and case["processes"]["thermodynamics"] == "dynamic":
        return 10.0
    return 600.0


def resolve(user_case: dict) -> dict:
    """Merge with defaults, fill derived defaults, validate; raises ValueError."""
    c = deep_merge(DEFAULTS, user_case)
    errs = []
    for path, allowed in ENUMS.items():
        try:
            v = get_path(c, path)
        except KeyError:
            continue
        if v not in allowed:
            errs.append(f"{path} must be one of {allowed} (got {v!r})")
    t = c["time"]
    if t["substep_s"] in (None, "", 0):
        t["substep_s"] = substep_default(c)
    for k in ("duration_s", "dt_s", "substep_s", "output_every_s"):
        if float(t[k]) <= 0 and not (k == "duration_s" and float(t[k]) == 0):
            errs.append(f"time.{k} must be > 0")
    if float(t["substep_s"]) > float(t["dt_s"]):
        t["substep_s"] = float(t["dt_s"])
    e = c["environment"]
    if not (150.0 < float(e["temperature_K"]) < 350.0):
        errs.append("environment.temperature_K out of range")
    if not (0.0 <= float(e["relative_humidity"]) <= 1.0):
        errs.append("environment.relative_humidity must be within [0, 1]")
    if float(e["pressure_Pa"]) <= 0:
        errs.append("environment.pressure_Pa must be > 0")
    if float((e.get("dilution") or {}).get("rate_s1", 0.0)) < 0:
        errs.append("environment.dilution.rate_s1 must be >= 0")
    gp = c["emissions"].get("gas_time_profile")
    if not (gp is None or gp == "jno2" or isinstance(gp, dict)):
        errs.append("emissions.gas_time_profile must be null, 'jno2' or {time_s, factor}")
    if c["kernel"]["isorropia"] == "none" and c["processes"]["wet_diameter"] == "isorropia":
        # no-ISORROPIA runs diagnose water with kappa-Koehler; keep the CESM setting
        pass
    if errs:
        raise ValueError("; ".join(errs))
    return c


def load(path: str | Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def save(case: dict, path: str | Path) -> None:
    Path(path).write_text(json.dumps(case, indent=2, ensure_ascii=False), encoding="utf-8")
