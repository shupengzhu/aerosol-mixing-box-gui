"""Generate the built-in scenario library (scenarios/*.json).

Run from the repository root:  python3 tools/make_scenarios.py
The scenarios are generated, not hand-written, so that every conversion from the
legacy SCRAM BoxApp (SCRAM 1.1 fork, H2O organics, 30 species) to the hetfix kernel
(14 species, CESM VBS organics) is explicit and reproducible.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "scenarios"

# --------------------------------------------------------------------------------------
# helpers
H = 3600.0


def lognormal_from_volume(dv_um, sigma, vol_um3_cm3):
    """Seigneur et al. (1986) modes are given as volume median diameter and volume;
    return number median diameter (um) and number (#/cm3) (Hatch-Choate)."""
    ls2 = math.log(sigma) ** 2
    dn = dv_um * math.exp(-3.0 * ls2)
    n = vol_um3_cm3 / (math.pi / 6.0 * dn ** 3 * math.exp(4.5 * ls2))
    return dn, n


HAZY = [(0.044, 1.2, 0.09), (0.24, 1.8, 5.8), (6.0, 2.2, 25.9)]     # Seigneur et al. 1986


def hazy_modes(fraction, comp, tag):
    out = []
    for name, (dv, s, v) in zip(("nuclei", "accumulation", "coarse"), HAZY):
        dn, n = lognormal_from_volume(dv, s, v * fraction)
        out.append({"name": f"{tag} {name}", "number_cm3": round(n, 6), "dg_um": round(dn, 6),
                    "sigma_g": s, "mass_fraction": comp})
    return out


def diurnal(hours, start_h, f):
    return [f((start_h + h) % 24.0) for h in hours]


def sun(h, peak):
    """Clear-sky-like daytime shape between 06 and 18 local time."""
    return peak * max(0.0, math.sin(math.pi * (h - 6.0) / 12.0)) ** 1.3


# --------------------------------------------------------------------------------------
# legacy organics -> CESM VBS (C* = 0.01, 0.1, 1, 10, 100 ug/m3 for SOA1P..SOA5P)
R, T0, TORR = 8.314, 298.0, 133.322
# partition coefficient Kp (m3/ug) or saturation pressure (torr) and MW (g/mol),
# from the legacy SCRAM ModuleInitialization.f90 (H2O module data)
LEGACY_ORG = {
    "BiA2D": {"kp": 6.25e-3, "mw": 186}, "BiA1D": {"kp": 2.73e-3, "mw": 170},
    "BiA0D": {"kp": 4.82e-5, "mw": 168}, "GLYOXAL": {"kp": 6.56e-4, "mw": 58},
    "MGLY": {"kp": 5.78e-12, "mw": 72}, "BiMT": {"kp": 0.8052, "mw": 136},
    "BiPER": {"kp": 0.1109, "mw": 168}, "BiDER": {"kp": 2.8, "mw": 136},
    "BiMGA": {"kp": 1.1281e-2, "mw": 120}, "AnBlP": {"p_torr": 6.8e-8, "mw": 167},
    "AnBmP": {"p_torr": 8.4e-6, "mw": 152}, "BiBlP": {"p_torr": 6.0e-10, "mw": 298},
    "BiBmP": {"p_torr": 3.0e-7, "mw": 236}, "BiNGA": {"p_torr": 1.39e-5, "mw": 165},
    "NIT3": {"p_torr": 1.45e-6, "mw": 272}, "BiNIT": {"p_torr": 2.5e-6, "mw": 215},
    "AnClP": {"kp": 55.56, "mw": 167}, "SOAlP": {"kp": 110.0, "mw": 392},
    "SOAmP": {"kp": 1.16, "mw": 392}, "SOAhP": {"kp": 0.031, "mw": 392},
    "POAlP": {"kp": 1.1, "mw": 280}, "POAmP": {"kp": 0.0116, "mw": 280},
    "POAhP": {"kp": 0.00031, "mw": 280},
}
CSTAR_MAX_GAS = 1000.0     # gas of species with C* above this is dropped (outside the VBS range)

# CESM-SCRAM sea-salt emission split (carma_model_mod, MAM4 ncl -> PNA/PHCL); Ca, K, Mg are
# diagnosed from PNA by the v2.3 kernel. CMAQ sea_spc (Na .3086, Cl .5538, SO4 .0776) is not
# used: without the 6 % Mg/Ca/K it is anion-rich and v1 outgasses ~25 % of the chloride.
SEA_SALT = {"PNA": 0.3866, "PHCL": 0.6134}


def cstar(name):
    d = LEGACY_ORG[name]
    if "kp" in d:
        return 1.0 / d["kp"]
    return d["p_torr"] * TORR * d["mw"] / (R * T0) * 1e6


def vbs_bin(name):
    k = round(math.log10(cstar(name))) + 3          # log10 C* = -2..2 -> 1..5
    return min(max(k, 1), 5)


# --------------------------------------------------------------------------------------
# legacy Greater Paris case (Zhu et al., 2015, GMD, section 4); core/templates/baseline12h.cfg
PARIS_BOUNDS = [0.001, 0.00316, 0.01, 0.0398, 0.1585, 0.6310, 2.5119, 10.0]
PARIS_NUMBER_M3 = [5895.412, 1138386.25, 5265435136.0, 8862283776.0, 492587168.0,
                   3598435.75, 13657.6796875]
# name: (gas_init ug/m3, gas_emission ug/m3/s, particle mass per bin ug/m3)
PARIS = {
    "MD": (0.0, 0.0, [5.91382257207e-16, 5.81481000503e-12, 6.5469424726e-07, 0.150394171476, 1.30960583687, 1.78676152229, 1.07002174854]),
    "BC": (0.0, 0.0, [3.70283163398e-15, 6.1401730922e-11, 0.0153744360432, 0.220083221793, 0.843515992165, 0.624725580215, 0.0312604494393]),
    "Na": (0.0, 0.0, [9.05962816105e-17, 4.45471296569e-13, 4.23234993718e-08, 0.00142551504541, 0.122185088694, 0.463657706976, 0.0248505901545]),
    "SO4": (0.000795035914052, 1.7216629566975E-05, [1.0621693412e-12, 2.51451020006e-08, 0.00195511663333, 0.165135219693, 3.2797396183, 1.08144080639, 0.0616511851549]),
    "NH4": (1.13932681084, 0.0, [8.60751921128e-12, 1.72809677679e-07, 0.0115803396329, 0.780958533287, 0.850924491882, 0.0343741178513, 0.00148804776836]),
    "NO3": (5.16737508774, 0.0, [4.08259970969e-12, 7.76404078806e-08, 0.00537681253627, 0.183640256524, 0.145689919591, 0.00693256594241, 0.00026601247373]),
    "Cl": (2.59052030256e-16, 0.0, [1.21632131281e-16, 7.96222408752e-13, 7.58435589887e-08, 0.0025532306172, 0.218745335937, 0.785844743252, 0.0420205518603]),
    "BiA2D": (0.19140727818, 0.0, [0.0] * 7),
    "BiA1D": (0.269935309887, 0.0, [1.70000420533e-13, 3.11097836381e-09, 0.000239889966906, 0.027067091316, 0.0356881357729, 0.00126388540957, 6.8415509304e-05]),
    "BiA0D": (0.000157813119586, 0.0, [2.95291308766e-12, 5.17785139209e-08, 0.00373264891095, 0.435509443283, 1.70179259777, 0.189896702766, 0.012012084946]),
    "GLYOXAL": (0.0, 0.0, [6.42989724752e-17, 5.11176654492e-16, 1.16171574568e-15, 2.13224145913e-14, 2.96884749857e-13, 8.90298983624e-14, 6.01279170158e-15]),
    "MGLY": (0.184095069766, 0.0, [6.43032407272e-17, 5.11187242404e-16, 7.71123040367e-15, 3.29417114246e-12, 6.42694925057e-12, 2.63978969977e-13, 1.44255911199e-14]),
    "BiMT": (1.0833561646e-05, 0.0, [9.02496532292e-17, 2.12309705407e-12, 3.27611957118e-07, 5.02633483848e-05, 0.00010233863577, 9.64058472164e-06, 5.80148878271e-07]),
    "BiPER": (0.105028316379, 0.0, [1.76526198173e-12, 3.59382426041e-08, 0.00230345316231, 0.123101621866, 0.113051958382, 0.00528930267319, 0.000262255896814]),
    "BiDER": (0.00204521627165, 0.0, [1.05224843053e-13, 3.96600530195e-09, 0.0004692550865, 0.0406863354146, 0.0682814195752, 0.00695510115474, 0.000396480987547]),
    "BiMGA": (0.0537621118128, 0.0, [5.97606680938e-14, 1.31010235993e-09, 8.93445612746e-05, 0.00474796351045, 0.00394762866199, 0.000126121929497, 5.06083006258e-06]),
    "AnBlP": (0.0278133861721, 0.0, [1.25574218382e-14, 3.83882980959e-10, 4.06167418987e-05, 0.0115467002615, 0.0190492141992, 0.0012389512267, 8.07932228781e-05]),
    "AnBmP": (0.38260102272, 0.0, [1.47543706853e-13, 1.64576430262e-09, 9.71314511844e-05, 0.00570566114038, 0.00536402780563, 0.000192818784853, 9.02282317838e-06]),
    "BiBlP": (5.97768703301e-05, 0.0, [2.04920089008e-14, 3.71545127997e-10, 2.72824399872e-05, 0.0034929565154, 0.0424189344049, 0.00832103379071, 0.000494618027005]),
    "BiBmP": (0.0430435314775, 0.0, [1.43726458686e-12, 1.3813746591e-08, 0.00097780150827, 0.0425603017211, 0.0346591919661, 0.00150242354721, 6.23350570095e-05]),
    "BiNGA": (0.0458301603794, 0.0, [5.32443428336e-15, 3.42310874091e-10, 2.74770427495e-05, 0.00589258922264, 0.00816686544567, 0.000400694290875, 2.47522584687e-05]),
    "NIT3": (0.107222780585, 0.0, [1.09209080208e-15, 1.55772103549e-10, 1.79746893991e-05, 0.00345024745911, 0.00481949746609, 0.000220006215386, 1.21009670693e-05]),
    "BiNIT": (1.25714254379, 0.0, [3.04545234517e-13, 6.72598376994e-09, 0.000576474645641, 0.100442916155, 0.182760760188, 0.0113450400531, 0.000771872350015]),
    "AnClP": (0.000179971626494, 0.0, [3.83301424522e-15, 8.27726151509e-11, 6.71521047479e-06, 0.00544697698206, 0.0192144121975, 0.00330412294716, 0.000208180616028]),
    "SOAlP": (0.000744042044971, 0.0, [3.64059436352e-14, 9.07149022567e-10, 8.32891528262e-05, 0.0934789478779, 0.499331772327, 0.11021065712, 0.00654964521527]),
    "SOAmP": (0.1370074898, 0.0, [5.13501621866e-13, 4.47616770316e-08, 0.00366603210568, 0.446798145771, 0.743557572365, 0.0539985448122, 0.00317199085839]),
    "SOAhP": (1.55686557293, 0.0, [4.00091028016e-12, 5.4051920273e-08, 0.00320525141433, 0.17241075635, 0.168752208352, 0.00730455154553, 0.00033628201345]),
    "POAlP": (0.15563082695, 0.00091374465264, [5.49477564143e-12, 1.07055612375e-07, 0.00780125614256, 0.556008994579, 0.852210462093, 0.0666050314903, 0.00415565958247]),
    "POAmP": (1.65659248829, 0.00116959311999, [5.69759734184e-13, 1.18432588181e-08, 0.000892996380571, 0.0659306943417, 0.0811726078391, 0.00356189976446, 0.000208040422876]),
    "POAhP": (2.42145466805, 0.00157164083794, [2.14457341181e-14, 4.44614206563e-10, 3.32924828399e-05, 0.00239930301905, 0.00291988137178, 0.000127840001369, 7.46147179598e-06]),
}
# particle emission rows (ug/m3/s per bin) apply to species 1 and 2 (MD, BC)
PARIS_PEMIS = {"MD": [0.0, 0.0, 0.0, 3.65498664905E-05, 7.3099732981125E-05, 7.3099732981125E-05, 0.0],
               "BC": [0.0, 0.0, 7.614434289275E-05, 0.00022843666957, 0.00030458284309, 0.00030458284309, 0.0]}
INORG = {"MD": ("PMD", None), "BC": ("PBC", None), "Na": ("PNA", None), "SO4": ("PSO4", "SULF"),
         "NH4": ("PNH4", "NH3"), "NO3": ("PNO3", "HNO3"), "Cl": ("PHCL", "HCL")}


def paris_case(name, coag, cond, nucl):
    mass, gas, gsrc, notes = {}, {}, {}, []

    def add(d, k, v):
        if isinstance(v, list):
            d[k] = [a + b for a, b in zip(d.get(k, [0.0] * len(v)), v)]
        else:
            d[k] = d.get(k, 0.0) + v

    vbs_map = {}
    for sp, (g0, ge, m) in PARIS.items():
        if sp in INORG:
            aer, gname = INORG[sp]
            add(mass, aer, m)
            if gname:
                add(gas, gname, g0)
                if ge:
                    add(gsrc, gname, ge)
            continue
        k = vbs_bin(sp)
        vbs_map[sp] = (round(cstar(sp), 4), k)
        if sp.startswith("POA"):
            add(mass, "PPOAlP", m)          # CESM POA is non-volatile
        else:
            add(mass, f"PSOA{k}P", m)
        if cstar(sp) <= CSTAR_MAX_GAS:
            add(gas, f"SOA{k}P", g0)
            if ge:
                add(gsrc, f"SOA{k}P", ge)
        elif g0 or ge:
            notes.append(f"{sp}: C*={cstar(sp):.3g} ug/m3 > {CSTAR_MAX_GAS:g}, gas phase "
                         f"({g0:.3g} ug/m3, emission {ge:.3g} ug/m3/s) not represented")
    pem = {"PMD": PARIS_PEMIS["MD"], "PBC": PARIS_PEMIS["BC"]}
    case = {
        "case": {"name": name, "reference": "Zhu, Sartelet and Seigneur (2015), GMD 8, 1595-1612, section 4"},
        "kernel": {"isorropia": "v23"},
        "environment": {"temperature_K": 279.553955078, "pressure_Pa": 79805.765625,
                        "relative_humidity": 0.745575904846},
        "time": {"duration_s": 12 * H, "dt_s": 600.0, "output_every_s": 1800.0},
        "mixing": {"preset": "NC6"},
        "size_bins": {"mode": "explicit", "bounds_um": PARIS_BOUNDS},
        "initial": {"particles": {"mode": "bin_mass", "mixing": "by_group_set",
                                  "sets": [["PMD"], ["PBC"],
                                           ["PNA", "PSO4", "PNH4", "PNO3", "PHCL"],
                                           ["PSOA1P", "PSOA2P", "PSOA3P", "PSOA4P", "PSOA5P", "PPOAlP"]],
                                  "number_cm3": [x / 1e6 for x in PARIS_NUMBER_M3],
                                  "mass_ugm3": {k: [round(x, 15) for x in v] for k, v in mass.items()
                                                if sum(v) > 0}},
                    "gas_ugm3": {k: v for k, v in gas.items() if v > 0}},
        "emissions": {"particles": {"mode": "bin_mass", "mixing": "by_group",
                                    "mass_ugm3_s": pem},
                      "gas_ugm3_s": gsrc},
        "processes": {"coagulation": coag, "condensation": cond,
                      "nucleation": {"enabled": nucl, "model": "ternary"}},
    }
    notes = ([f"Legacy H2O organic species mapped to the CESM VBS bin nearest to log10 C* "
              f"(298 K): " + ", ".join(f"{s}->{'PPOAlP' if s.startswith('POA') else f'PSOA{k}P'} "
                                       f"(C*={c:g})" for s, (c, k) in vbs_map.items()),
              "Particulate POA is mapped to the non-volatile CESM POA (PPOAlP); POA gas and its "
              "emissions are mapped to the VBS gas bin of the same C* (SVOC).",
              "Initial particles: legacy SCRAM external start; dust, BC, inorganic and organic "
              "populations are separate classes, number split by volume.",
              "Legacy solver ROS2 and redistribution 'moving diameter' replaced by the CESM-SCRAM "
              "defaults (ETR, moving-center dual pivot)."] + notes)
    return case, notes


# --------------------------------------------------------------------------------------
def scenario(sid, category, name_en, name_zh, sum_en, sum_zh, case, notes=(), tags=()):
    case = dict(case)
    case.setdefault("case", {})
    case["case"].update({"name": name_en, "scenario_id": sid, "category": category})
    return {"id": sid, "category": category, "name": {"en": name_en, "zh": name_zh},
            "summary": {"en": sum_en, "zh": sum_zh}, "tags": list(tags), "notes": list(notes),
            "case": case}


def build():
    out = []
    # ---- teaching
    out.append(scenario(
        "tutorial_minimal", "teaching",
        "Minimal tutorial (BC + sulfate coagulation)", "最简教学案例（BC + 硫酸盐凝并）",
        "Fresh BC and aged sulfate populations coagulate for one hour; watch number move from "
        "the pure classes into mixed classes.",
        "新鲜 BC 与老化硫酸盐两个颗粒群凝并 1 小时，观察数浓度从纯组分类群进入混合类群。",
        {"kernel": {"isorropia": "none"},
         "environment": {"temperature_K": 298.0, "relative_humidity": 0.5},
         "time": {"duration_s": 1 * H, "dt_s": 300.0, "output_every_s": 300.0},
         "mixing": {"preset": "NC6"},
         "size_bins": {"mode": "log", "n": 8, "dmin_um": 0.01, "dmax_um": 2.5},
         "initial": {"particles": {"mode": "lognormal", "modes": [
             {"name": "fresh BC", "number_cm3": 20000, "dg_um": 0.04, "sigma_g": 1.6,
              "mass_fraction": {"PBC": 1.0}},
             {"name": "sulfate", "number_cm3": 5000, "dg_um": 0.12, "sigma_g": 1.5,
              "mass_fraction": {"PSO4": 1.0}}]}},
         "processes": {"coagulation": True, "condensation": False}},
        tags=["coagulation", "BC", "mixing state"]))

    # ---- GMD validation (Zhu et al. 2015, section 3): hazy distribution, half sulfate half BC
    hazy = hazy_modes(0.5, {"PSO4": 1.0}, "sulfate") + hazy_modes(0.5, {"PBC": 1.0}, "BC")
    for sid, coag, en, zh in (("gmd_hazy_condensation", False, "GMD hazy condensation validation",
                               "GMD hazy 冷凝验证"),
                              ("gmd_hazy_coag_cond", True, "GMD hazy coagulation + condensation validation",
                               "GMD hazy 凝并+冷凝验证")):
        out.append(scenario(
            sid, "gmd_validation", en, zh,
            "Hazy tri-modal distribution (Seigneur et al., 1986), externally mixed half sulfate / "
            "half BC by volume, sulfuric acid produced at 2.29e-4 ug/m3/s for 12 h, 298 K, 1 atm, "
            "RH 70%" + (", with coagulation." if coag else "."),
            "Seigneur et al. (1986) hazy 三模态分布，按体积一半硫酸盐一半 BC 外混，硫酸气体以 "
            "2.29e-4 ug/m3/s 生成 12 小时，298 K、1 atm、RH 70%" + ("，含凝并。" if coag else "。"),
            {"case": {"reference": "Zhu, Sartelet and Seigneur (2015), GMD 8, 1595-1612, section 3"},
             "kernel": {"isorropia": "none"},
             "environment": {"temperature_K": 298.0, "pressure_Pa": 101325.0, "relative_humidity": 0.7},
             "time": {"duration_s": 12 * H, "dt_s": 600.0, "output_every_s": 1800.0},
             "mixing": {"preset": "NC6"},
             "size_bins": {"mode": "log", "n": 4, "dmin_um": 0.01, "dmax_um": 10.0},
             "initial": {"particles": {"mode": "lognormal", "modes": hazy}},
             "emissions": {"gas_ugm3_s": {"SULF": 2.29e-4}},
             "processes": {"coagulation": coag, "condensation": True}},
            notes=["Legacy template used 4 size sections and 2 composition sections; here 4 log "
                   "sections 0.01-10 um and the NC6 classes (sulfate and BC fall in separate classes). "
                   "Increase the number of sections for a converged reference.",
                   "Hazy modes converted from volume median diameter and volume to number median "
                   "diameter and number (Hatch-Choate)."],
            tags=["validation", "condensation", "sulfate", "BC"]))

    # ---- GMD Greater Paris A-D
    for sid, letter, coag, cond, nucl, en_p, zh_p in (
            ("gmd_paris_emission_only", "A", False, False, False, "emission only", "仅排放"),
            ("gmd_paris_coagulation", "B", True, False, False, "emission + coagulation", "排放 + 凝并"),
            ("gmd_paris_condensation", "C", False, True, False, "emission + condensation", "排放 + 冷凝"),
            ("gmd_paris_full", "D", True, True, True, "full dynamics", "全过程")):
        case, notes = paris_case(f"GMD Greater Paris scenario {letter}", coag, cond, nucl)
        out.append(scenario(
            sid, "gmd_reference", f"GMD Greater Paris scenario {letter} ({en_p})",
            f"GMD 巴黎场景 {letter}（{zh_p}）",
            f"Greater Paris reference case of Zhu et al. (2015), scenario {letter}: {en_p}. "
            "Polyphemus-derived initial particles in 7 sections, BC and dust emissions, sulfuric "
            "acid and organic vapour sources, 12 h.",
            f"Zhu et al. (2015) 巴黎参考场景 {letter}：{zh_p}。Polyphemus 提取的 7 个粒径段初始颗粒、"
            "BC 与沙尘排放、硫酸与有机气体源，12 小时。",
            case, notes=notes, tags=["Paris", "emission", "reference"]))

    # ---- mixing-state scenarios (docs/case_format.md section 5.1)
    hours = [float(h) for h in range(0, 25)]
    start = 6.0
    T = diurnal(hours, start, lambda h: 293.0 + 6.0 * math.cos(math.pi * (h - 14.0) / 12.0))
    RH = diurnal(hours, start, lambda h: 0.70 - 0.22 * math.cos(math.pi * (h - 14.0) / 12.0))
    J = diurnal(hours, start, lambda h: sun(h, 8.0e-3))
    emis = [1.0 if h < 12.0 else 0.0 for h in hours]
    out.append(scenario(
        "urban_plume_aging", "mixing_state", "Urban plume aging", "城市羽流老化",
        "Lagrangian urban air parcel from 06:00 for 24 h (after the PartMC urban-plume case): "
        "background Aitken and accumulation modes, diesel, gasoline and cooking emissions for "
        "12 h, dilution with background air, photochemical production of sulfuric acid, nitric "
        "acid and SOA vapours following J_NO2. BC-containing particles age from fresh to coated.",
        "仿照 PartMC urban plume 的拉格朗日城市气团，06:00 起 24 小时：背景 Aitken 与积聚模态，"
        "柴油车、汽油车与烹饪排放持续 12 小时，与背景空气稀释混合，硫酸、硝酸与 SOA 前体按 J_NO2 "
        "光化学生成。含 BC 颗粒由新鲜逐步被包覆老化。",
        {"case": {"reference": "Riemer, West, Zaveri and Easter (2009), JGR 114, D09202; "
                               "PartMC scenarios/urban_plume"},
         "kernel": {"isorropia": "v23"},
         "environment": {"temperature_K": T[0], "pressure_Pa": 101325.0, "relative_humidity": RH[0],
                         "time_series": {"time_s": [h * H for h in hours], "temperature_K": T,
                                         "relative_humidity": RH, "jno2_s1": J},
                         "dilution": {"rate_s1": 1.5e-5, "background": "initial"}},
         "time": {"duration_s": 24 * H, "dt_s": 600.0, "output_every_s": 1800.0},
         "mixing": {"preset": "NC6"},
         "size_bins": {"mode": "log", "n": 16, "dmin_um": 0.005, "dmax_um": 2.5},
         "initial": {"particles": {"mode": "lognormal", "modes": [
             {"name": "background Aitken", "number_cm3": 3200, "dg_um": 0.02, "sigma_g": 1.45,
              "mass_fraction": {"PSO4": 0.364, "PNH4": 0.136, "PSOA3P": 0.5}},
             {"name": "background accumulation", "number_cm3": 2900, "dg_um": 0.116, "sigma_g": 1.65,
              "mass_fraction": {"PSO4": 0.364, "PNH4": 0.136, "PSOA3P": 0.5}}]},
             "gas_ugm3": {"NH3": 5.0, "HNO3": 1.0, "SOA3P": 0.5, "SOA4P": 1.0, "SOA5P": 2.0}},
         "emissions": {"particles": {"mode": "lognormal", "modes": [
             {"name": "diesel", "number_cm3_s": 0.4, "dg_um": 0.05, "sigma_g": 1.7,
              "mass_fraction": {"PBC": 0.7, "PPOAlP": 0.3}},
             {"name": "gasoline", "number_cm3_s": 0.125, "dg_um": 0.05, "sigma_g": 1.7,
              "mass_fraction": {"PBC": 0.2, "PPOAlP": 0.8}},
             {"name": "meat cooking", "number_cm3_s": 0.0225, "dg_um": 0.086, "sigma_g": 1.9,
              "mass_fraction": {"PPOAlP": 1.0}}]},
             "gas_ugm3_s": {"SULF": 2.5e-4, "HNO3": 5.0e-4, "SOA1P": 1.5e-5, "SOA2P": 3.0e-5,
                            "SOA3P": 4.5e-5, "SOA4P": 9.0e-5, "SOA5P": 1.2e-4},
             "time_profile": {"time_s": [h * H for h in hours], "factor": emis},
             "gas_time_profile": "jno2"}},
        notes=["Particle emission fluxes of the PartMC urban plume (m-2 s-1) divided by a 400 m "
               "mixing height; emissions stop after 12 h.",
               "SCRAM has no gas-phase chemistry: the condensable vapours are prescribed sources "
               "scaled by J_NO2/max(J_NO2) (noon rates SULF 0.9, HNO3 1.8, SOA vapours 1.1 ug/m3/h).",
               "Dilution 1.5e-5 1/s towards the initial (background) state."],
        tags=["urban", "BC aging", "mixing state", "dilution", "photochemistry"]))

    out.append(scenario(
        "seasalt_nitrate", "mixing_state", "Sea salt + nitrate (chloride depletion)",
        "海盐 + 硝酸盐（氯亏损）",
        "Coastal air: fresh sea spray (CMAQ sea-salt speciation) meets polluted continental air "
        "with ammonium sulfate particles and nitric acid. HNO3 displaces chloride from sea salt "
        "(NaCl + HNO3 -> NaNO3 + HCl); ISORROPIA II also sees the sea-salt Ca, K, Mg.",
        "沿海气团：新鲜海盐（CMAQ 海盐组分）与含硫酸铵颗粒和硝酸的污染大陆气团混合。HNO3 置换海盐中的 "
        "氯（NaCl + HNO3 -> NaNO3 + HCl）；ISORROPIA II 同时考虑海盐中的 Ca、K、Mg。",
        {"case": {"reference": "CESM-SCRAM sea-salt split; CMAQ v5.5 AERO_DATA sea_spc cation ratios"},
         "kernel": {"isorropia": "v23"},
         "environment": {"temperature_K": 293.0, "pressure_Pa": 101325.0, "relative_humidity": 0.80,
                         "time_series": {"time_s": [h * H for h in hours[:13]],
                                         "jno2_s1": [sun(10.0 + h, 7.0e-3) for h in hours[:13]]}},
         "time": {"duration_s": 12 * H, "dt_s": 600.0, "output_every_s": 1800.0},
         "mixing": {"preset": "NC6"},
         "size_bins": {"mode": "log", "n": 14, "dmin_um": 0.01, "dmax_um": 10.0},
         "initial": {"particles": {"mode": "lognormal", "modes": [
             {"name": "sea salt film/jet", "number_cm3": 30, "dg_um": 0.25, "sigma_g": 2.0,
              "mass_fraction": SEA_SALT},
             {"name": "sea salt coarse", "number_cm3": 3, "dg_um": 1.5, "sigma_g": 2.0,
              "mass_fraction": SEA_SALT},
             {"name": "continental sulfate", "number_cm3": 1500, "dg_um": 0.08, "sigma_g": 1.6,
              "mass_fraction": {"PSO4": 0.727, "PNH4": 0.273}}]},
             "gas_ugm3": {"HNO3": 4.0, "NH3": 0.5}},
         "emissions": {"gas_ugm3_s": {"HNO3": 2.0e-4}, "gas_time_profile": "jno2"}},
        notes=["Sea-salt composition as emitted by CESM-SCRAM (PNA 0.3866, PHCL 0.6134). Ca, K and Mg "
               "are not tracers: the v2.3 kernel derives them from Na with the CMAQ ratios "
               "(kernel/data/speciation_cmaq_aero.json); v1 and no-ISORROPIA ignore them.",
               "Run with kernel v1 and v23 to see the effect of the crustal cations."],
        tags=["sea salt", "nitrate", "chloride depletion", "ISORROPIA II"]))

    rh_ramp = [0.30 + 0.55 * min(h / 10.0, 1.0) for h in hours[:13]]
    t_ramp = [300.0 - 0.6 * h for h in hours[:13]]
    out.append(scenario(
        "dust_crustal", "mixing_state", "Dust + crustal ions", "沙尘 + 地壳离子",
        "A dust plume (CMAQ dust speciation: Ca, K, Mg, Na as fractions of dust) mixes with "
        "polluted air rich in nitric acid and ammonia while RH rises from 30% to 85%. With "
        "ISORROPIA II the crustal cations neutralise nitrate on dust (Ca(NO3)2); with v1 dust is inert.",
        "沙尘气团（CMAQ 沙尘组分：Ca、K、Mg、Na 占沙尘质量比例）与富含硝酸和氨的污染气团混合，RH 从 30% "
        "升至 85%。ISORROPIA II 中地壳阳离子在沙尘上中和硝酸盐（Ca(NO3)2）；v1 中沙尘为惰性。",
        {"case": {"reference": "CMAQ v5.5 AERO_DATA dust_spc / asoil factors; Fountoukis and Nenes (2007) ACP"},
         "kernel": {"isorropia": "v23"},
         "environment": {"temperature_K": 300.0, "pressure_Pa": 101325.0, "relative_humidity": 0.30,
                         "time_series": {"time_s": [h * H for h in hours[:13]],
                                         "temperature_K": t_ramp, "relative_humidity": rh_ramp}},
         "time": {"duration_s": 12 * H, "dt_s": 600.0, "output_every_s": 1800.0},
         "mixing": {"preset": "NC6"},
         "size_bins": {"mode": "log", "n": 14, "dmin_um": 0.01, "dmax_um": 10.0},
         "initial": {"particles": {"mode": "lognormal", "modes": [
             {"name": "dust fine", "number_cm3": 20, "dg_um": 0.4, "sigma_g": 1.8,
              "mass_fraction": {"PMD": 1.0}},
             {"name": "dust coarse", "number_cm3": 12, "dg_um": 1.2, "sigma_g": 1.8,
              "mass_fraction": {"PMD": 1.0}},
             {"name": "polluted sulfate", "number_cm3": 2000, "dg_um": 0.08, "sigma_g": 1.6,
              "mass_fraction": {"PSO4": 0.727, "PNH4": 0.273}}]},
             "gas_ugm3": {"HNO3": 6.0, "NH3": 4.0, "HCL": 0.5}},
         "emissions": {"gas_ugm3_s": {"SULF": 1.0e-4}}},
        notes=["Ca, K, Mg, Na are diagnosed by the v2.3 kernel from PMD with the CMAQ fine (<2.5 um) "
               "and coarse fractions (kernel/data/speciation_cmaq_aero.json).",
               "Run with kernel v1 and v23: the difference is the crustal-ion effect."],
        tags=["dust", "crustal ions", "nitrate", "ISORROPIA II"]))
    return out


def main():
    OUT.mkdir(exist_ok=True)
    for s in build():
        (OUT / f"{s['id']}.json").write_text(json.dumps(s, indent=1, ensure_ascii=False) + "\n",
                                             encoding="utf-8")
        print("wrote", s["id"])


if __name__ == "__main__":
    main()
