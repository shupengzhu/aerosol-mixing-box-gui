"""Cost-accuracy analysis of SCRAM configurations against the particle-resolved reference.

python opt_analysis.py <screen_dir> <runs_dir> <strata_summary.json> <regimes_summary.json> <strata15.json> <out_dir>

screen_dir : runs/opt/screen/<id>.json (representation, PartMC projected on each configuration)
runs_dir   : runs/opt/runs/<id>.json (first batch) and runs/opt/runs/<id>/<cfg>.json
Writes: opt_long.csv (one row per state x configuration), opt_configs.csv (statistics per
configuration and subset), opt_stats.md.
"""
import glob, json, os, sys
from pathlib import Path

import numpy as np
import pandas as pd

SCR, RUNS, SS, RS, ST, OUT = sys.argv[1:7]
os.makedirs(OUT, exist_ok=True)
SS_ = ["0.1", "0.2", "0.5"]; WS = ["0.2", "0.5", "1.0"]
NSPEC = 14                                        # aerosol species incl. water (tracer count = sections x (species + number))

REG_TYPE = {"north_china_plain_east_china": "polluted land", "indo_gangetic_plain": "polluted land",
            "southeast_asia": "polluted land", "western_europe": "polluted land", "eastern_us": "polluted land",
            "amazon": "biomass burning", "central_africa_biomass": "biomass burning",
            "sahara_sahel_dust": "dust", "east_asian_dust": "dust", "middle_east_dust": "dust",
            "nw_pacific_outflow": "ocean", "north_atlantic": "ocean", "southern_ocean": "ocean",
            "southeast_atlantic_stratocumulus": "ocean", "arctic": "remote"}

# ---------------- weights and metadata ----------------
st = json.load(open(ST))
wst = {(s["region"], s["season"]): s for s in st["strata"]}
meta = {}
for s in st["samples"]:
    sid = f"strat_{s['region']}_{s['season']}_{s['kind']}"
    w = wst[(s["region"], s["season"])]
    meta[sid] = {"layer": "strata", "region": s["region"], "season": s["season"], "kind": s["kind"],
                 "type": REG_TYPE[s["region"]], "wa": w["area_time_frac"] / 2, "wc": w["lowcloud_frac"] / 2}
reg = json.load(open(RS))
for r in reg:
    meta[r["id"]] = {"layer": "regime", "region": "regime", "season": "-", "kind": "typical",
                     "type": "global background", "wa": r.get("weight"), "wc": r.get("weight_lowcloud")}
summ = {r["id"]: r for r in json.load(open(SS))}
summ.update({r["id"]: r for r in reg})

# ---------------- representation (screen) ----------------
rows = {}
for f in glob.glob(f"{SCR}/*.json"):
    d = json.load(open(f)); sid = d["id"]
    if sid.endswith("_internal") or sid not in meta:
        continue
    tr = d["truth"]
    for cfg, m in d["configs"].items():
        if "error" in m:
            continue
        r = rows.setdefault((sid, cfg), {"id": sid, "cfg": cfg})
        r["nsize"], r["ncls"] = m["nsize"], m["nclass"]
        for s in SS_:
            r[f"ccn_t{s}"] = tr["ccn"][s]; r[f"R{s}"] = m["ccn"][s] / tr["ccn"][s]
        for w in WS:
            r[f"nd_t{w}"] = tr["nd"][w]; r[f"Rd{w}"] = m["act"][w]["parcel"] / tr["nd"][w]
        if "eabs" in m:
            r["eabs_rep"], r["kbt_rep"], r["kce_rep"] = m["eabs"], m["k_bt"], m["k_cesm"]

# ---------------- process runs ----------------
def run_items():
    for f in glob.glob(f"{RUNS}/*.json"):
        d = json.load(open(f))
        for cfg, m in d["configs"].items():
            yield d["id"], cfg, m
    for f in glob.glob(f"{RUNS}/*/*.json"):
        d = json.load(open(f))
        for cfg, m in d["configs"].items():
            yield d["id"], cfg, m


for sid, cfg, m in run_items():
    if "error" in m or sid not in meta:
        continue
    r = rows.setdefault((sid, cfg), {"id": sid, "cfg": cfg})
    r["nsize"], r["ncls"] = m["nsize"], m["nclass"]
    r["cpu_s"] = m.get("cpu_s")
    tr = summ.get(sid, {}).get("truth", {})
    for s in SS_:
        t = r.get(f"ccn_t{s}") or tr.get("ccn", {}).get(s)
        if t:
            r[f"S{s}"] = m["ccn"][s] / t
        r[f"TS{s}"] = (m.get("time_mean_abs") or {}).get(s)
    for w in WS:
        t = r.get(f"nd_t{w}") or ((summ.get(sid, {}).get("activation") or {}).get(w, {}).get("truth") or {}).get("nd")
        a = m["act"][w]
        if t:
            r[f"Sd{w}"] = a["parcel"] / t                         # aerosol error in N_d (parcel on simulated state)
            for k in ("arg", "fn", "mbn"):
                if k in a:
                    r[f"E_{k}{w}"] = a[k] / t                  # end to end: scheme on simulated state / parcel on particles
                    r[f"A_{k}{w}"] = a[k] / a["parcel"]        # activation factor
            r[f"Dpos{w}"] = a.get("fn_D_pos")
    het = summ.get(sid, {}).get("het", {}).get("truth", {})
    et = tr.get("eabs")
    r["eabs_sim"] = m.get("eabs"); r["eabs_t"] = et
    if het:
        r["kbt_sim"] = m.get("k_bt") / het["bertram_thornton"] if m.get("k_bt") else None
        r["kce_sim"] = m.get("k_cesm") / het["cesm"] if m.get("k_cesm") else None
    r["bc"] = (summ.get(sid, {}).get("sim_NC6") or {}).get("bc")

# production runs reused (no N_d chain): NC1, NC11 at U1 and NC6 x8, NC9BC x1
for sid, s in summ.items():
    if sid not in meta:
        continue
    tr = s.get("truth", {})
    for lab, cfg, nc in (("sim_NC1", "NC1@U1", 1), ("sim_NC11", "NC11@U1", 11)):
        if lab in s:
            r = rows.setdefault((sid, cfg), {"id": sid, "cfg": cfg, "nsize": 16, "ncls": nc})
            for q in SS_:
                r.setdefault(f"S{q}", s[lab]["ccn"][q] / tr["ccn"][q])
            if "cpu_s" in s[lab]:
                r.setdefault("cpu_prod", s[lab].get("cpu_s"))
    for lab, cfg, ns, nc in (("NC6_x8", "NC6@U8", 128, 6), ("NC9BC_x1", "NC9BC@U1", 16, 9)):
        nm = (s.get("numerics") or {}).get(lab)
        if nm:
            r = rows.setdefault((sid, cfg), {"id": sid, "cfg": cfg, "nsize": ns, "ncls": nc})
            for q in SS_:
                r.setdefault(f"S{q}", nm["scram"]["ccn"][q] / tr["ccn"][q])
            r.setdefault("eabs_sim", nm["scram"].get("eabs_prog")); r.setdefault("eabs_t", tr.get("eabs"))
            r.setdefault("cpu_prod", nm["scram"].get("cpu_s"))

L = pd.DataFrame(rows.values())
M = pd.DataFrame.from_dict(meta, orient="index")
L = L.join(M, on="id")
L["scheme"] = L.cfg.str.split("@").str[0]; L["grid"] = L.cfg.str.split("@").str[1]
L["tracers"] = L.nsize * L.ncls * (NSPEC + 1)
for s in SS_:
    if f"S{s}" in L and f"R{s}" in L:
        L[f"P{s}"] = L[f"S{s}"] / L[f"R{s}"]
# composite: mean over the three supersaturations of |S - 1| (end of run) and of the 24 h time-mean |error|
L["Sall"] = 1 + L[[f"S{s}" for s in SS_]].sub(1).abs().mean(axis=1, skipna=False)
L["TSall"] = 1 + L[[f"TS{s}" for s in SS_]].mean(axis=1, skipna=False) if all(f"TS{s}" in L for s in SS_) else np.nan
L["Eall_mbn"] = 1 + L[[f"E_mbn{w}" for w in WS]].sub(1).abs().mean(axis=1, skipna=False) if all(f"E_mbn{w}" in L for w in WS) else np.nan
L["Sdall"] = 1 + L[[f"Sd{w}" for w in WS]].sub(1).abs().mean(axis=1, skipna=False) if all(f"Sd{w}" in L for w in WS) else np.nan
# relative box cost: median over states of cpu(cfg)/cpu(NC6@U1) for the same state
if len(sys.argv) > 7:                             # per-host normalised costs (cpu / cpu of NC6@U1, same host, same state)
    cr = pd.read_csv(sys.argv[7]).set_index(["id", "cfg"])["cost_rel"]
    L["cost_rel"] = [cr.get((i, c), np.nan) for i, c in zip(L.id, L.cfg)]
else:
    base = L[L.cfg == "NC6@U1"].set_index("id")["cpu_s"]
    L["cost_rel"] = L.cpu_s / L.id.map(base)


L.to_csv(f"{OUT}/opt_long.csv", index=False)


def wq(x, w, q):
    x, w = np.asarray(x, float), np.asarray(w, float)
    ok = np.isfinite(x) & np.isfinite(w)
    x, w = x[ok], w[ok]
    if len(x) == 0:
        return np.nan
    o = np.argsort(x); x, w = x[o], w[o]
    c = np.cumsum(w) / w.sum()
    return float(np.interp(q, c, x))


def stats(g, col, wcol=None):
    x = (g[col] - 1).abs() * 100
    if wcol:
        return wq(x, g[wcol], 0.5), wq(x, g[wcol], 0.9)
    x = x.dropna()
    return (float(x.median()), float(x.quantile(0.9))) if len(x) else (np.nan, np.nan)


SUBSETS = {"strata": lambda d: d[d.layer == "strata"],
           "strata_lowcloud": lambda d: d[d.layer == "strata"],
           "polluted_land": lambda d: d[d.type == "polluted land"],
           "regimes_area": lambda d: d[d.layer == "regime"]}
WCOL = {"strata": None, "strata_lowcloud": "wc", "polluted_land": None, "regimes_area": "wa"}
cols = [f"R{s}" for s in SS_] + [f"S{s}" for s in SS_] + [f"P{s}" for s in SS_] + [f"Rd{w}" for w in WS] + \
       [f"Sd{w}" for w in WS] + [f"E_{k}{w}" for k in ("arg", "fn", "mbn") for w in WS] + ["kbt_sim", "kce_sim"] + \
       ["Sall", "TSall", "Eall_mbn", "Sdall"]
out = []
for cfg, g0 in L.groupby("cfg"):
    for sub, f in SUBSETS.items():
        g = f(g0)
        if g.empty:
            continue
        row = {"cfg": cfg, "subset": sub, "n": len(g), "ncls": g.ncls.iloc[0], "nsize": g.nsize.iloc[0],
               "tracers": g.tracers.iloc[0], "cost_rel": float(g.cost_rel.median()) if g.cost_rel.notna().any() else np.nan}
        for c in cols:
            if c in g and g[c].notna().sum() >= 5:
                row[c + "_med"], row[c + "_p90"] = stats(g, c, WCOL[sub])
        out.append(row)
C = pd.DataFrame(out)
C.to_csv(f"{OUT}/opt_configs.csv", index=False)
print(C[C.subset == "strata"].sort_values("tracers")[["cfg", "n", "ncls", "nsize", "tracers", "cost_rel"] +
      [c for c in ("R0.1_p90", "S0.1_med", "S0.1_p90", "S0.2_p90", "Sd0.5_p90", "E_mbn0.5_med", "E_mbn0.5_p90") if c in C]].round(2).to_string())
