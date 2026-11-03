"""Particle-resolved and MAM4-modal N2O5 uptake rate (Bertram-Thornton 2009 and CESM fixed gamma) at
the end of a PartMC run, with exactly the scrambox hetchem formulas: every particle is one population.
Usage: python tools/het_truth.py <scenario_run_dir> [...]   -> <dir>/het_truth.json"""
import json, math, sys
from pathlib import Path
import numpy as np
T_ = Path(__file__).resolve().parent; ROOT = T_.parent; sys.path.insert(0, str(T_)); sys.path.insert(0, str(ROOT))
from scrambox import case as casemod
from scrambox.compile import build_grid
from scrambox.partmc import read_particles, _pmc_name
from scrambox.post import hetchem
from scrambox.post.common import Populations
from wp2_mam4_projection import assign_modes, particles as mam_particles
OPTS = {"gases": ["N2O5"], "rules": ["cesm", "bertram_thornton"]}
def kvals(res):
    g = res["rules"]["N2O5"] if "rules" in res else res["N2O5"]
    out = {}
    for rule in ("cesm", "bertram_thornton"):
        v = g["resolved"][rule] if "resolved" in g else g[rule]
        out[rule] = v.get("k_s1"); out[rule + "_gamma_eff"] = v.get("gamma_eff")
    return out
for d in map(Path, sys.argv[1:]):
    pdir = d / "partmc_s1"
    c = casemod.load(pdir / "case.resolved.json"); g = build_grid(c); table = g["table"]
    sp = list(table.aerosol); inv = {_pmc_name(s): i for i, s in enumerate(sp)}
    dens = np.array([table.density_gcm3(s) for s in sp])
    last = None
    for rec in read_particles(pdir):
        last = rec
    t, env, spec, mass, nconc, gnames, gppb = last
    m = np.zeros((len(sp), mass.shape[1]))
    for k, s in enumerate(spec):
        if s in inv: m[inv[s]] = mass[k]
    keep = (m.sum(0) - m[sp.index("PH2O")]) > 0
    num = nconc[keep][:, None]; ms = (m[:, keep] * nconc[keep] * 1e9).T[:, None, :]          # ug m-3 per particle population
    T, P = env.get("temperature", 298.0), env.get("pressure", 1e5)
    pop = Populations(number=num, mass=ms, species=sp, density_gcm3=dens, bounds_um=np.array([0.001, 100.0]),
                      class_labels=["particle"], T=T, RH=0.0, P=P)
    res = {"truth": kvals(hetchem.compute(pop, OPTS))}
    # MAM4 static modal representation: species summed within each MAM4 mode (internal mixture)
    rec = None
    for rec in mam_particles(pdir):
        pass
    _, T2, P2, names, vol, n, kap = rec
    mode = assign_modes(names, vol, kap)
    from wp2_mam4_projection import mam4_metrics
    from wp2_activation_parcel import lognormal_sections
    if len(mode) == int(keep.sum()):
        mpops = mam4_metrics(T2, P2, names, vol, n, kap, mode, True)["modes"]
        NUMS, MASSES = [], []
        for p in mpops:
            sel = mode == p["mode"]; mtot = ms[sel, 0, :].sum(0)                 # ug m-3 of each species in the mode
            dry = mtot.sum() - mtot[sp.index("PH2O")]
            if dry <= 0: continue
            nn, rd = lognormal_sections(p["N_cm3"] * 1e6, p["dgn_um"], p["sigma"], nsec=60)   # MAM4 lognormal, fixed sigma
            vfrac = nn * rd ** 3 / (nn * rd ** 3).sum()                        # dry-volume share of each sub-section
            for j in range(len(nn)):
                if nn[j] <= 0: continue
                NUMS.append([nn[j]]); MASSES.append(mtot * vfrac[j])
        popm = Populations(number=np.array(NUMS), mass=np.array(MASSES)[:, None, :], species=sp, density_gcm3=dens,
                           bounds_um=np.array([0.001, 100.0]), class_labels=["mode"], T=T, RH=0.0, P=P)
        res["mam4"] = kvals(hetchem.compute(popm, OPTS)); res["mam4_modes"] = [p["mode"] for p in mpops]
    else:
        res["mam4_error"] = f"particle count mismatch {len(mode)} vs {int(keep.sum())}"
    (d / "het_truth.json").write_text(json.dumps(res, indent=1, default=float)); print(d.name, res, flush=True)
