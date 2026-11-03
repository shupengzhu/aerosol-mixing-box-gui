"""SCRAM species tables: read the kernel's species files and write new ones.

The kernel reads its species table from SCRAM_SPECIES_FILE (scram_species_io_mod).
The library files in kernel/data are the SCRAM runtime tables:
  - scram_species_iso16_20260731.dat  ISORROPIA profile (14 aerosol species, VBS organics)
  - scram_species_noiso_hetfix_nc6.dat  no-ISORROPIA profile (11 aerosol species)
A case may add inert custom species; they are written into the generated table.
"""
from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "kernel" / "data"
LIBRARY_FILES = {
    "iso": DATA_DIR / "scram_species_iso16_20260731.dat",
    "noiso": DATA_DIR / "scram_species_noiso_hetfix_nc6.dat",
}

# Sections holding "<name> <value>" pairs (possibly several pairs per line).
PAIR_SECTIONS = [
    "mass_density_aer", "scram_id_aer", "molecular_weight_aer", "molecular_diameter_aer",
    "collision_factor_aer", "saturation_pressure", "partition_coefficient",
    "vaporization_enthalpy", "accomodation_coefficient", "surface_tension",
]
# Sections holding one "<name> <name>" relation per line.
RELATION_SECTIONS = ["aerosol_species_group_relations", "gas_species_aerosol_interact"]
# Sections holding a plain list of names.
LIST_SECTIONS = ["species", "aerosol_species", "soluble_species_aer", "aerosol_groups",
                 "isorropia_species", "pankow_species", "poa_species"]
SECTION_ORDER = ["species", "aerosol_species", "soluble_species_aer", "aerosol_groups",
                 "aerosol_species_group_relations", "mass_density_aer", "scram_id_aer",
                 "isorropia_species", "pankow_species", "poa_species",
                 "gas_species_aerosol_interact", "molecular_weight_aer",
                 "molecular_diameter_aer", "collision_factor_aer", "saturation_pressure",
                 "partition_coefficient", "vaporization_enthalpy",
                 "accomodation_coefficient", "surface_tension"]


@dataclass
class SpeciesTable:
    """In-memory form of a SCRAM species file."""
    gases: list[str]
    aerosol: list[str]                      # in scram id order, PH2O last
    soluble: list[str]
    groups: list[str]
    group_of: dict[str, str]                # aerosol species -> group (PH2O has none)
    isorropia: list[str]
    pankow: list[str]
    poa: list[str]
    gas_of: dict[str, str]                  # aerosol species -> condensing gas
    props: dict[str, dict[str, float]] = field(default_factory=dict)   # section -> {species: value}

    # ---- derived helpers
    def density_gcm3(self, name: str) -> float:
        return float(self.props["mass_density_aer"][name])

    def mw_gmol(self, name: str) -> float:
        return float(self.props["molecular_weight_aer"].get(name, 100.0))

    def group_index(self, name: str) -> int:
        """1-based kernel group index; 0 for water (excluded from dry groups)."""
        if name == "PH2O":
            return 0
        return self.groups.index(self.group_of[name]) + 1

    def gas_index(self, name: str) -> int:
        return self.gases.index(name)

    def to_dict(self) -> dict:
        return {
            "gases": self.gases, "aerosol": self.aerosol, "groups": self.groups,
            "group_of": self.group_of, "isorropia": self.isorropia, "pankow": self.pankow,
            "poa": self.poa, "gas_of": self.gas_of,
            "density_gcm3": {s: self.props["mass_density_aer"].get(s) for s in self.aerosol},
            "mw_gmol": {s: self.props["molecular_weight_aer"].get(s) for s in self.aerosol},
        }


def _tokens(line: str) -> list[str]:
    return line.split("#", 1)[0].split()


def read_species_file(path: str | Path) -> SpeciesTable:
    sections: dict[str, list[list[str]]] = {}
    current = None
    for raw in Path(path).read_text().splitlines():
        line = raw.strip()
        m = re.match(r"^\[(\w+)\]$", line)
        if m:
            current = m.group(1)
            sections.setdefault(current, [])
            continue
        toks = _tokens(line)
        if current and toks:
            sections[current].append(toks)

    def flat(sec):
        return [t for row in sections.get(sec, []) for t in row]

    props: dict[str, dict[str, float]] = {}
    for sec in PAIR_SECTIONS:
        vals = flat(sec)
        props[sec] = {vals[i]: float(vals[i + 1]) for i in range(0, len(vals) - 1, 2)}
    rel = {row[0]: row[1] for row in sections.get("aerosol_species_group_relations", []) if len(row) >= 2}
    gas_of = {row[0]: row[1] for row in sections.get("gas_species_aerosol_interact", []) if len(row) >= 2}
    aerosol = flat("aerosol_species")
    ids = props["scram_id_aer"]
    aerosol = sorted(aerosol, key=lambda s: ids.get(s, 1e9))
    return SpeciesTable(gases=flat("species"), aerosol=aerosol, soluble=flat("soluble_species_aer"),
                        groups=flat("aerosol_groups"), group_of=rel, isorropia=flat("isorropia_species"),
                        pankow=flat("pankow_species"), poa=flat("poa_species"), gas_of=gas_of,
                        props=props)


def library(profile: str) -> SpeciesTable:
    """profile: 'iso' (ISORROPIA v1 or v2.3) or 'noiso'."""
    return read_species_file(LIBRARY_FILES[profile])


def with_custom_species(table: SpeciesTable, custom: list[dict]) -> SpeciesTable:
    """Add inert custom aerosol species (inserted before PH2O, which stays last).
    Each item: {"name", "group", "density_gcm3", "mw_gmol"}; group must exist."""
    t = copy.deepcopy(table)
    for c in custom:
        name = c["name"].strip()
        if not re.match(r"^[A-Za-z][A-Za-z0-9_]{0,15}$", name):
            raise ValueError(f"invalid species name '{name}'")
        if name in t.aerosol:
            raise ValueError(f"species '{name}' already exists")
        if c["group"] not in t.groups:
            raise ValueError(f"unknown group '{c['group']}' for species '{name}'")
        water_at = t.aerosol.index("PH2O")
        t.aerosol.insert(water_at, name)
        t.soluble.append(name)
        t.group_of[name] = c["group"]
        t.props["mass_density_aer"][name] = float(c["density_gcm3"])
        t.props["molecular_weight_aer"][name] = float(c.get("mw_gmol", 100.0))
    for i, s in enumerate(t.aerosol, start=1):
        t.props["scram_id_aer"][s] = float(i)
    return t


def _pairs(names, values, per_line=5, fmt=lambda v: repr(v)):
    items = [f"{n} {fmt(values[n])}" for n in names if n in values]
    return [" ".join(items[i:i + per_line]) for i in range(0, len(items), per_line)]


def _num(v: float) -> str:
    if float(v).is_integer() and abs(v) < 1e6:
        return f"{int(v)}" if abs(v) >= 1 else f"{v:g}"
    return f"{v:.6g}"


def write_species_file(t: SpeciesTable, path: str | Path) -> None:
    out = []
    for sec in SECTION_ORDER:
        out.append(f"[{sec}]")
        if sec == "species":
            out += [" ".join(t.gases[i:i + 10]) for i in range(0, len(t.gases), 10)]
        elif sec == "aerosol_species":
            out.append(" ".join(t.aerosol))
        elif sec == "soluble_species_aer":
            out.append(" ".join(s for s in t.aerosol if s in t.soluble))
        elif sec == "aerosol_groups":
            out.append(" ".join(t.groups))
        elif sec == "aerosol_species_group_relations":
            out += [f"{s} {t.group_of[s]}" for s in t.aerosol if s in t.group_of]
        elif sec == "isorropia_species":
            out.append(" ".join(t.isorropia))
        elif sec == "pankow_species":
            out.append(" ".join(t.pankow))
        elif sec == "poa_species":
            out.append(" ".join(t.poa))
        elif sec == "gas_species_aerosol_interact":
            out += [f"{s} {t.gas_of[s]}" for s in t.aerosol if s in t.gas_of]
        elif sec == "scram_id_aer":
            out += _pairs(t.aerosol, t.props[sec], per_line=14, fmt=lambda v: str(int(v)))
        else:
            out += _pairs(t.aerosol, t.props[sec], fmt=_num)
        out.append("")
    Path(path).write_text("\n".join(out))
