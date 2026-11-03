"""Letter figures, revision 3 (ensemble-based).
Usage: python make_figs_v3.py <ensemble_table.csv> <results_json dir> <dust_diag.json> <out_dir>
Fig 1  error budget: R, P (and MAM4) for CCN over the 126+30 scenarios by region type; k_N2O5 and E_abs panels.
Fig 2  mechanisms: (a) size-refinement convergence for the three mechanism cases and the nitrate-growth case decomposition,
       (b) cost-accuracy of size vs class refinement (urban), (c) nitrate-growth deficit: variants and resolution convergence.
Fig 3  activation: parameterized/parcel N_d for ARG and FN across the ensemble, split by the FN branch; error chain.
"""
import sys, json, os
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
CSV, RJ, DD, OUT = sys.argv[1:5]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6, "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                     "legend.frameon": False, "savefig.dpi": 300, "axes.spines.top": False, "axes.spines.right": False, "xtick.direction": "out", "ytick.direction": "out"})
C = {"mam4": "#9e9e9e", "rep": "#3b7bb5", "sim": "#d9622b", "x8": "#f3b48a", "arg": "#8c8c8c", "fn": "#1f7a5a", "fnD": "#7fbf9f", "ref": "#222222"}
T = pd.read_csv(CSV); ext = T[(T.layer == "strata") & (~T.internal)]; reg = T[T.layer == "regime"]; allr = pd.concat([ext, reg])
order = ["polluted land", "biomass burning", "dust", "ocean", "remote", "global background"]
labels = {"polluted land": "polluted land\n(n=%d)", "biomass burning": "biomass\nburning (n=%d)", "dust": "dust (n=%d)", "ocean": "ocean (n=%d)", "remote": "remote (n=%d)", "global background": "global\nbackground (n=%d)"}
pct = lambda x: 100 * (np.asarray(x, float) - 1)
def letter(ax, t, dx=-0.12):
    ax.text(dx, 1.03, t, transform=ax.transAxes, fontweight="bold", fontsize=9, ha="left", va="bottom")
rng = np.random.default_rng(3)
def strip(ax, cols, ylab=True, xlim=(-60, 250), band=5):
    for ti, typ in enumerate(order):
        sel = allr[allr.type == typ]
        for j, (k, col, lab) in enumerate(cols):
            v = pct(sel[k]); v = v[np.isfinite(v)]
            if len(v) == 0: continue
            y = ti * (len(cols) + 1) + j
            ax.scatter(v, y + rng.uniform(-0.28, 0.28, len(v)), s=5, color=col, alpha=0.55, lw=0, label=lab if ti == 0 else None, zorder=2)
            ax.plot([np.median(v)] * 2, [y - 0.45, y + 0.45], color="k", lw=1.3, zorder=3)
    ax.axvline(0, color="k", lw=0.6); ax.axvspan(-band, band, color="0.92", zorder=0, lw=0)
    ax.set_yticks([ti * (len(cols) + 1) + (len(cols) - 1) / 2 for ti in range(len(order))])
    ax.set_yticklabels([labels[t] % len(allr[allr.type == t]) for t in order] if ylab else [""] * len(order))
    ax.invert_yaxis(); ax.set_xscale("symlog", linthresh=10, linscale=1.0); ax.set_xlim(*xlim)
    ax.set_xticks([-50, -10, 0, 10, 50, 200]); ax.set_xticklabels(["-50", "-10", "0", "10", "50", "200"]); ax.tick_params(axis="y", length=0)
# ---------------------------------------------------------------- Figure 1
fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), gridspec_kw={"width_ratios": [1.25, 1, 1]})
cols = [("M0.1", C["mam4"], "modal representation (MAM4 modes)"), ("R0.1", C["rep"], "sectional representation (six classes)"),
        ("S0.1", C["sim"], "sectional simulation"), ("X8_0.1", C["x8"], "simulation, 8× size sections")]
strip(axs[0], cols); axs[0].set_title("CCN at 0.1% supersaturation", fontsize=8, loc="left")
colsk = [("K_mam", C["mam4"], None), ("K_rep", C["rep"], None), ("K_sim", C["sim"], None)]
strip(axs[1], colsk, ylab=False); axs[1].set_title("N$_2$O$_5$ uptake rate $k$ (BT09)", fontsize=8, loc="left")
E = allr.copy(); E = E[(E.bc > 0.1) & (~E.type.isin(["dust", "ocean"]))]
E["ER"] = E.eabs_rep / E.eabs_t; E["ES"] = E.eabs_sim / E.eabs_t
save = allr; allr = E
strip(axs[2], [("ER", C["rep"], None), ("ES", C["sim"], None)], ylab=False, xlim=(-30, 250)); allr = save
axs[2].set_title("BC absorption enhancement", fontsize=8, loc="left")
axs[1].set_xlabel("Bias relative to particle-resolved reference (%)")
h, l = axs[0].get_legend_handles_labels(); fig.legend(h, l, loc="upper center", ncol=4, bbox_to_anchor=(0.5, 1.0), markerscale=2.4, handletextpad=0.2, columnspacing=1.0)
for ax, t in zip(axs, "abc"): letter(ax, t, dx=-0.55 if t == "a" else -0.1)
fig.tight_layout(rect=(0, 0, 1, 0.92), w_pad=0.8); fig.savefig(f"{OUT}/fig1_error_budget.png"); fig.savefig(f"{OUT}/fig1_error_budget.pdf"); plt.close(fig)
# ---------------------------------------------------------------- Figure 2
e1 = json.load(open(f"{RJ}/e1.json")); e2 = json.load(open(f"{RJ}/e2.json")); d3 = json.load(open(f"{RJ}/diag3.json")); dd = json.load(open(DD)); e5 = json.load(open(f"{RJ}/e5.json"))
fig, axs = plt.subplots(1, 3, figsize=(7.2, 2.6))
ax = axs[0]  # cost vs accuracy, urban
tr = e1["truth"]
for mix, col, mk, lab in (("NC6", C["sim"], "o", "six classes"), ("NC9BC", "#6a51a3", "D", "nine BC-resolved classes")):
    xs, yc, ye = [], [], []
    for k in (1, 2, 4, 8):
        s = e1[f"scram_{mix}_x{k}"]; xs.append(s["cpu_s"] / e1["scram_NC6_x1"]["cpu_s"]); yc.append(100 * abs(s["ccn0.1"] / tr["ccn0.1"] - 1)); ye.append(100 * abs(s["eabs"] / 1.77 - 1))
    ax.plot(xs, yc, "-", marker=mk, color=col, ms=4, lw=1, label=f"{lab}: CCN(0.1%)")
    ax.plot(xs, ye, "--", marker=mk, color=col, ms=4, lw=1, mfc="white", label=f"{lab}: $E_{{\\rm abs}}$")
    if mix == "NC6":
        for x, k, y in zip(xs, (1, 2, 4, 8), yc): ax.annotate(f"×{k}", (x, y), xytext=(4, 3), textcoords="offset points", fontsize=6)
ax.set_xscale("log"); ax.set_xlabel("Cost relative to six classes, base grid"); ax.set_ylabel("|bias| (%)"); ax.set_ylim(0, 38); ax.set_xlim(0.7, 1000)
ax.legend(fontsize=6, loc="upper center", handlelength=2.2, bbox_to_anchor=(0.62, 1.0)); ax.set_title("Urban plume, 24 h", fontsize=7, loc="left")
ax.set_ylim(0, 48)
ax = axs[1]  # nitrate-growth deficit: convergence and variants
conv = [(1, d3["moving_center_dualpivot_x1"]["ccn0.1"]), (8, d3["moving_center_dualpivot_x8"]["ccn0.1"]), (16, d3["moving_center_dualpivot_x16"]["ccn0.1"]), (32, d3["moving_center_dualpivot_x32"]["ccn0.1"])]
tru = dd["full"]["truth"]["0.1"][-1]
ax.plot([4.7 * k for k, _ in conv], [v / 1487.1 for _, v in conv], "o-", color=C["sim"], ms=4, lw=1, label="dust + crustal ions (full)")
for name, col, lab in (("no_hno3", "#1f7a5a", "no HNO$_3$"), ("no_dust", "#b8860b", "no dust"), ("nocoag", "#6a51a3", "no coagulation")):
    v = dd[name]; ax.plot([4.7, 37.6], [v["scram_x1"]["0.1"][-1] / v["truth"]["0.1"][-1], v["scram_x8"]["0.1"][-1] / v["truth"]["0.1"][-1]], "s--", color=col, ms=3.5, lw=1, label=lab)
ax.axhline(1, color=C["ref"], lw=0.7, ls="--"); ax.set_xscale("log", base=2); ax.set_xticks([5, 10, 20, 40, 80, 150]); ax.set_xticklabels(["5", "10", "20", "40", "80", "150"])
ax.set_xlabel("Size sections per decade"); ax.set_ylabel("CCN(0.1%) / reference"); ax.set_ylim(0.55, 1.08); ax.legend(fontsize=6, loc="lower right")
ax.set_title("Nitrate growth case, one class", fontsize=7, loc="left")
ax = axs[2]  # integrators, six classes
for solver, col in (("euler", "#6a51a3"), ("etr", C["sim"]), ("ros2", C["rep"])):
    pts = sorted((float(k.split("_sub")[1]), v["ccn0.1"] / 1487.1) for k, v in e2.items() if f"NC6_{solver}_sub" in k and "default" not in k)
    ax.plot([p[0] for p in pts], [p[1] for p in pts], "o-", color=col, ms=3.5, lw=1, label={"euler": "Euler", "etr": "ETR", "ros2": "ROS2"}[solver])
ax.axhline(1, color=C["ref"], lw=0.7, ls="--"); ax.set_xscale("log"); ax.set_ylim(0.8, 1.06); ax.set_xlabel("Integrator sub-step (s)"); ax.set_ylabel("CCN(0.1%) / reference")
ax.legend(fontsize=6.5, loc="lower right"); ax.set_title("Six classes, 8× sections", fontsize=7, loc="left")
for ax, t in zip(axs, "abc"): letter(ax, t, dx=-0.3)
fig.tight_layout(w_pad=1.0); fig.savefig(f"{OUT}/fig2_mechanisms.png"); fig.savefig(f"{OUT}/fig2_mechanisms.pdf"); plt.close(fig)
# ---------------------------------------------------------------- Figure 3
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [1.15, 1]})
ax = axs[0]
rows = []
for w in ("0.2", "0.5", "1.0"):
    for _, r in ext.iterrows():
        rows.append({"w": w, "arg": r[f"A_arg{w}"], "fn": r[f"A_fn{w}"], "D": r[f"A_fnD{w}"], "rep": r[f"A_rep{w}"], "proc": r[f"A_proc{w}"]})
A = pd.DataFrame(rows)
xpos = {"0.2": 0, "0.5": 1, "1.0": 2}
for w in ("0.2", "0.5", "1.0"):
    a = A[A.w == w]; x0 = xpos[w] * 4
    for j, (k, col, lab, sel) in enumerate((("rep", C["rep"], "representation", slice(None)), ("proc", C["sim"], "process", slice(None)), ("arg", C["arg"], "ARG", slice(None)), ("fn", C["fn"], "FN, analytic branch", a.D == True), ("fn", C["fnD"], "FN, empirical branch", a.D == False))):
        v = a[k][sel] if not isinstance(sel, slice) else a[k]; v = v[np.isfinite(v)]
        x = x0 + (j - 2) * 0.55
        ax.scatter(np.full(len(v), x) + rng.uniform(-0.15, 0.15, len(v)), v, s=5, color=col, alpha=0.5, lw=0, label=lab if w == "0.2" else None, zorder=2)
        ax.plot([x - 0.22, x + 0.22], [np.median(v)] * 2, color="k", lw=1.3, zorder=3)
ax.axhline(1, color="k", lw=0.6); ax.axhspan(0.9, 1.1, color="0.92", zorder=0, lw=0); ax.set_yscale("log"); ax.set_ylim(0.1, 2.5)
ax.set_yticks([0.125, 0.25, 0.5, 1, 2]); ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2"]); ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax.set_xticks([0, 4, 8]); ax.set_xticklabels(["0.2", "0.5", "1"]); ax.set_xlabel("Updraft $w$ (m s$^{-1}$)")
ax.set_ylabel("Factor applied to $N_d$"); ax.legend(fontsize=6, loc="lower left", ncol=2, markerscale=2.2, handletextpad=0.1, columnspacing=0.8)
ax.set_title("Droplet-number error chain, %d scenarios" % len(ext), fontsize=7, loc="left")
ax = axs[1]  # FN vs ARG total factor vs reference parcel N_d, colored by branch
tot_fn = []; tot_arg = []
for w in ("0.2", "0.5", "1.0"):
    for _, r in ext.iterrows():
        f = r[f"A_rep{w}"] * r[f"A_proc{w}"]; tot_fn.append((f * r[f"A_fn{w}"], r[f"A_fnD{w}"])); tot_arg.append(f * r[f"A_arg{w}"])
tf = np.array([t[0] for t in tot_fn]); dpos = np.array([t[1] == True for t in tot_fn]); ta = np.array(tot_arg)
bins = np.geomspace(0.1, 2.5, 28)
ax.hist(ta[np.isfinite(ta)], bins=bins, color=C["arg"], alpha=0.6, label=f"ARG (median {np.nanmedian(ta):.2f})")
ax.hist(tf[np.isfinite(tf) & dpos], bins=bins, color=C["fn"], alpha=0.7, label=f"FN, analytic (median {np.nanmedian(tf[dpos]):.2f})")
ax.hist(tf[np.isfinite(tf) & ~dpos], bins=bins, color=C["fnD"], alpha=0.8, label=f"FN, empirical (median {np.nanmedian(tf[~dpos]):.2f})")
ax.axvline(1, color="k", lw=0.6); ax.set_xscale("log"); ax.set_xticks([0.125, 0.25, 0.5, 1, 2]); ax.set_xticklabels(["0.125", "0.25", "0.5", "1", "2"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xlabel("$N_d$ on simulated state / $N_d$ on particles (parcel)"); ax.set_ylabel("Scenario × updraft count"); ax.legend(fontsize=6, loc="upper left")
ax.set_title("End-to-end droplet-number factor", fontsize=7, loc="left")
for ax, t in zip(axs, "ab"): letter(ax, t, dx=-0.16)
fig.tight_layout(w_pad=1.2); fig.savefig(f"{OUT}/fig3_activation.png"); fig.savefig(f"{OUT}/fig3_activation.pdf"); plt.close(fig)
json.dump({"fig3_end_to_end": {"arg_median": float(np.nanmedian(ta)), "fn_analytic_median": float(np.nanmedian(tf[dpos])), "fn_empirical_median": float(np.nanmedian(tf[~dpos])), "n_D_neg": int((~dpos).sum()), "n": int(len(tf))}},
          open(f"{OUT}/fig3_data.json", "w"), indent=1)
print("done")
