"""Ensemble analysis for the letter: error-budget factors R (representation), P (process), A (activation)
over the stratified (112 + 14 internal) and area-regime (30) scenarios.
Usage: python ensemble_figs.py strata_summary.json regimes_summary.json strata.json coast110.json out_dir"""
import json, sys, os, re
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
S, RG, STJ, COAST, OUT = sys.argv[1:6]
os.makedirs(OUT, exist_ok=True)
strata = json.load(open(S)); regimes = json.load(open(RG)); st = json.load(open(STJ)); coast = json.load(open(COAST))
plt.rcParams.update({"font.size": 8, "axes.linewidth": 0.6, "legend.frameon": False, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False})
C = {"mam4": "#b0b0b0", "rep": "#3b7bb5", "sim": "#d9622b", "x8": "#f0a070", "arg": "#8c8c8c", "fn": "#1f7a5a"}
REG_TYPE = {"north_china_plain_east_china": "polluted land", "indo_gangetic_plain": "polluted land", "southeast_asia": "polluted land",
            "western_europe": "polluted land", "eastern_us": "polluted land", "amazon": "biomass burning", "central_africa_biomass": "biomass burning",
            "sahara_sahel_dust": "dust", "east_asian_dust": "dust", "middle_east_dust": "dust", "nw_pacific_outflow": "ocean", "north_atlantic": "ocean",
            "southern_ocean": "ocean", "southeast_atlantic_stratocumulus": "ocean", "arctic": "remote"}
# stratum weights
W = {(x["region"], x["season"]): x for x in st["strata"]}
def ratio(a, b, floor=1.0):
    return a / b if (a is not None and b not in (None, 0) and b > floor) else np.nan
rows = []
for r in strata + regimes:
    if "truth" not in r or "error" in r: continue
    sid = r["id"]; internal = sid.endswith("_internal")
    if sid.startswith("strat_"):
        reg, sea, kind = r["region"], r["season"], r["kind"]; w = W[(reg, sea)]
        wa, wp, wc = w["area_time_frac"] / 2, w["pop_time_frac"] / 2, w["lowcloud_frac"] / 2; typ = REG_TYPE[reg]; layer = "strata"
    else:
        reg, sea, kind = "regime", "-", "typical"; wa, wc = r["weight"], r["weight_lowcloud"]; wp = r.get("pop_weight", np.nan); typ = "global background"; layer = "regime"
    t = r["truth"]
    row = {"id": sid, "layer": layer, "internal": internal, "region": reg, "season": sea, "kind": kind, "type": typ, "wa": wa, "wp": wp, "wc": wc,
           "N": t["N"], "lat": r["lat"], "lon": r["lon"], "pm": r.get("PM_fine_ugm3")}
    for s in ("0.1", "0.2", "0.5"):
        tc = t["ccn"][s]
        row[f"R{s}"] = ratio(r["rep_NC6"]["ccn"][s], tc); row[f"P{s}"] = ratio(r["sim_NC6"]["ccn"][s], r["rep_NC6"]["ccn"][s]) if tc > 1 else np.nan
        row[f"S{s}"] = ratio(r["sim_NC6"]["ccn"][s], tc); row[f"M{s}"] = ratio(r.get("mam4", {}).get("ccn", {}).get(s), tc)
        row[f"R1_{s}"] = ratio(r["rep_NC1"]["ccn"][s], tc); row[f"S1_{s}"] = ratio(r["sim_NC1"]["ccn"][s], tc)
        n = r.get("numerics", {})
        row[f"X8_{s}"] = ratio(n.get("NC6_x8", {}).get("scram", {}).get("ccn", {}).get(s), tc)
        row[f"B9_{s}"] = ratio(n.get("NC9BC_x1", {}).get("scram", {}).get("ccn", {}).get(s), tc)
        row[f"TR{s}"] = r["rep_NC6"]["time"][s]["mean_abs"]; row[f"TS{s}"] = r["sim_NC6"]["time"][s]["mean_abs"]
    row["chi_t"] = t["chi"]; row["chi_rep"] = r["rep_NC6"]["chi"]; row["chi_sim"] = r["sim_NC6"]["chi"]
    row["eabs_t"] = t.get("eabs"); row["eabs_rep"] = r["rep_NC6"].get("eabs"); row["eabs_sim"] = r["sim_NC6"].get("eabs"); row["bc"] = r["sim_NC6"].get("bc")
    h = r.get("het") or {}
    if h and h.get("truth"):
        kt = h["truth"]["bertram_thornton"]; kc = h["truth"]["cesm"]
        row["K_rep"] = ratio(r["rep_NC6"].get("k_bt"), kt, 1e-7); row["K_sim"] = ratio(r["sim_NC6"].get("k_bt"), kt, 1e-7)
        row["K_mam"] = ratio((h.get("mam4") or {}).get("bertram_thornton"), kt, 1e-7)
        row["K_cesm_rep"] = ratio(r["rep_NC6"].get("k_cesm"), kc, 1e-7); row["K_cesm_sim"] = ratio(r["sim_NC6"].get("k_cesm"), kc, 1e-7)
        row["K_cesm_mam"] = ratio((h.get("mam4") or {}).get("cesm"), kc, 1e-7)
    a = r.get("activation", {})
    for w_ in ("0.2", "0.5", "1.0"):
        if w_ in a:
            x = a[w_]; tr = x["truth"]["nd"]
            row[f"A_rep{w_}"] = ratio(x["p_nc6"]["parcel"], tr); row[f"A_proc{w_}"] = ratio(x["scram_nc6"]["parcel"], x["p_nc6"]["parcel"]) if tr > 1 else np.nan
            row[f"A_arg{w_}"] = ratio(x["scram_nc6"]["arg"], x["scram_nc6"]["parcel"]); row[f"A_fn{w_}"] = ratio(x["scram_nc6"]["fn"], x["scram_nc6"]["parcel"])
            row[f"A_fnD{w_}"] = x["scram_nc6"].get("fn_D_pos"); row[f"A_mam{w_}"] = ratio(x["mam4_parcel"]["nd"], tr)
    rows.append(row)
import csv
keys = sorted({k for r in rows for k in r}, key=lambda k: (k not in ("id", "layer", "region", "season", "kind", "type"), k))
with open(f"{OUT}/ensemble_table.csv", "w", newline="") as f:
    wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader(); [wr.writerow(r) for r in rows]
def wq(x, w, q):
    x = np.asarray(x, float); w = np.asarray(w, float); m = np.isfinite(x) & np.isfinite(w) & (w > 0)
    if m.sum() == 0: return np.nan
    o = np.argsort(x[m]); xs, ws = x[m][o], w[m][o]; cw = np.cumsum(ws) / ws.sum(); return float(np.interp(q, cw, xs))
def pct(x): return 100 * (np.asarray(x, float) - 1)
ext = [r for r in rows if r["layer"] == "strata" and not r["internal"]]
reg_rows = [r for r in rows if r["layer"] == "regime"]
# ---------------- summary statistics table (markdown) --------------------------------------------
L = ["# Ensemble error budget (CCN, end of run). Entries: median |bias| %, weighted P90 |bias| % (weights: area, population, low cloud; unweighted for the stratified sample)", ""]
def block(sel, name, wkey=None):
    L.append(f"## {name} (n={len(sel)})"); L.append("| quantity | R (representation) | P (process, sim/rep) | S = R*P (simulation) | MAM4 modal | NC6 x8 | NC9BC |"); L.append("|---|---|---|---|---|---|---|")
    for s in ("0.1", "0.2", "0.5"):
        cells = []
        for k in (f"R{s}", f"P{s}", f"S{s}", f"M{s}", f"X8_{s}", f"B9_{s}"):
            v = np.abs(pct([r[k] for r in sel])); w = np.array([r[wkey] for r in sel]) if wkey else np.ones(len(sel))
            cells.append(f"{wq(v, w, 0.5):.1f} / {wq(v, w, 0.9):.1f}" if np.isfinite(v).any() else "-")
        L.append(f"| CCN({s}%) | " + " | ".join(cells) + " |")
    for lab, k in (("k_N2O5 (BT09) rep", "K_rep"), ("k_N2O5 sim", "K_sim"), ("k_N2O5 MAM4", "K_mam")):
        v = np.abs(pct([r.get(k, np.nan) for r in sel])); w = np.array([r[wkey] for r in sel]) if wkey else np.ones(len(sel))
        if np.isfinite(v).any(): L.append(f"| {lab} | {wq(v, w, 0.5):.1f} / {wq(v, w, 0.9):.1f} | | | | | |")
    L.append("")
block(ext, "Stratified sample, unweighted")
block(ext, "Stratified sample, population-weighted", "wp")
block(ext, "Stratified sample, low-cloud-weighted", "wc")
block(reg_rows, "Area regimes (global background), area-weighted", "wa")
for typ in sorted({r["type"] for r in ext}):
    block([r for r in ext if r["type"] == typ], f"Stratified: {typ}")
block([r for r in ext if r["kind"] == "polluted"], "Stratified: polluted days only")
L.append("## Activation factor A (parameterized N_d / parcel N_d on the SCRAM state), median / P10-P90, stratified sample")
L.append("| w (m/s) | ARG | FN | FN where analytic branch (D>=0) | FN where empirical (D<0) | n D<0 | representation factor | process factor |"); L.append("|---|---|---|---|---|---|---|---|")
for w_ in ("0.2", "0.5", "1.0"):
    arg = np.array([r.get(f"A_arg{w_}", np.nan) for r in ext]); fn = np.array([r.get(f"A_fn{w_}", np.nan) for r in ext]); D = np.array([r.get(f"A_fnD{w_}") for r in ext], dtype=object)
    dp = np.array([d is True for d in D], bool); dn = np.array([d is False for d in D], bool)
    if len(ext) == 0: continue
    f = lambda x: f"{np.nanmedian(x):.2f} ({np.nanpercentile(x, 10):.2f}-{np.nanpercentile(x, 90):.2f})" if (len(x) and np.isfinite(np.asarray(x, float)).any()) else "-"
    L.append(f"| {w_} | {f(arg)} | {f(fn)} | {f(fn[dp])} | {f(fn[dn])} | {dn.sum()} | {f(np.array([r.get(f'A_rep{w_}', np.nan) for r in ext]))} | {f(np.array([r.get(f'A_proc{w_}', np.nan) for r in ext]))} |")
L.append("")
# internal vs external init
ints = {r["id"].replace("_internal", ""): r for r in rows if r["internal"]}
if ints:
    L.append("## Internal-mixture initialisation sensitivity (14 polluted days): simulation bias S and MAM4 bias, external vs internal init"); L.append("| id | S0.1 ext | S0.1 int | M0.1 ext | M0.1 int | R0.1 ext | R0.1 int |"); L.append("|---|---|---|---|---|---|---|")
    for k, ri in ints.items():
        re_ = next((r for r in ext if r["id"] == k), None)
        if re_: L.append(f"| {k.replace('strat_', '')} | {pct(re_['S0.1']):+.1f} | {pct(ri['S0.1']):+.1f} | {pct(re_['M0.1']):+.1f} | {pct(ri['M0.1']):+.1f} | {pct(re_['R0.1']):+.1f} | {pct(ri['R0.1']):+.1f} |")
open(f"{OUT}/ensemble_stats.md", "w", encoding="utf-8").write("\n".join(L))
# ---------------- Figure 1: error-budget distributions --------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6), sharey=False)
order = ["polluted land", "biomass burning", "dust", "ocean", "remote", "global background"]
allrows = ext + reg_rows
for ax, (s, title) in zip(axs, (("0.1", "CCN at 0.1% supersaturation"), ("0.2", "CCN at 0.2%"), ("0.5", "CCN at 0.5%"))):
    y0 = 0
    for ti, typ in enumerate(order):
        sel = [r for r in allrows if r["type"] == typ]
        if not sel: continue
        for j, (k, col, lab) in enumerate(((f"M{s}", C["mam4"], "modal (MAM4) representation"), (f"R{s}", C["rep"], "sectional representation"), (f"S{s}", C["sim"], "sectional simulation"), (f"X8_{s}", C["x8"], "simulation, 8x size sections"))):
            v = pct([r[k] for r in sel]); v = v[np.isfinite(v)]
            if len(v) == 0: continue
            y = ti * 5 + j
            ax.scatter(v, np.full(len(v), y) + np.random.default_rng(ti * 10 + j).uniform(-0.3, 0.3, len(v)), s=5, color=col, alpha=0.6, lw=0, label=lab if ti == 0 else None)
            ax.plot([np.median(v)] * 2, [y - 0.45, y + 0.45], color="k", lw=1.2)
    ax.axvline(0, color="k", lw=0.6); ax.axvspan(-5, 5, color="0.93", zorder=0, lw=0)
    ax.set_yticks([ti * 5 + 1.5 for ti in range(len(order))]); ax.set_yticklabels(order if ax is axs[0] else [""] * len(order))
    ax.set_title(title, fontsize=8, loc="left"); ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=10); ax.set_xticks([-50, -10, 0, 10, 50, 200]); ax.set_xticklabels(["-50", "-10", "0", "10", "50", "200"])
axs[1].set_xlabel("Bias relative to particle-resolved reference (%)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=4, fontsize=7, bbox_to_anchor=(0.5, 1.0), markerscale=2.5)
fig.tight_layout(rect=(0, 0, 1, 0.9)); fig.savefig(f"{OUT}/figE1_error_budget.png"); fig.savefig(f"{OUT}/figE1_error_budget.pdf"); plt.close(fig)
# ---------------- Figure S1: sample map --------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.2, 3.4))
for sg in coast:
    a = np.array(sg); ax.plot(a[:, 0], a[:, 1], color="0.5", lw=0.4)
cols = {"polluted land": "#c0392b", "biomass burning": "#e67e22", "dust": "#b8860b", "ocean": "#3b7bb5", "remote": "#7f8c8d", "global background": "#2ecc71"}
for typ in order:
    sel = [r for r in allrows if r["type"] == typ]
    x = [(r["lon"] - 360 if r["lon"] > 180 else r["lon"]) for r in sel]; y = [r["lat"] for r in sel]
    m = ["^" if r["kind"] == "polluted" else "o" for r in sel]
    for xi, yi, mi in zip(x, y, m): ax.plot(xi, yi, mi, color=cols[typ], ms=4 if mi == "o" else 5, mfc=cols[typ] if mi == "^" else "white", mew=0.8)
    ax.plot([], [], "o", color=cols[typ], label=f"{typ} (n={len(sel)})")
ax.plot([], [], "^", color="k", label="polluted day (P90)"); ax.plot([], [], "o", color="k", mfc="white", label="typical day")
ax.set_xlim(-180, 180); ax.set_ylim(-90, 90); ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.legend(fontsize=6, loc="upper center", ncol=4, bbox_to_anchor=(0.5, -0.12), frameon=False)
fig.tight_layout(); fig.savefig(f"{OUT}/figS1_sample_map.png"); plt.close(fig)
print("\n".join(L[:40]))
