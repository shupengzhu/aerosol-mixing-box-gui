"""Fig. 4 (cost-accuracy) of letter v4.
python opt_figs.py <opt_configs.csv> <opt_long.csv> <out_dir>
(a) CCN simulation error (90th percentile over the stratified states) against box cost,
(b) the same, low-cloud weighted, (c) droplet number end to end (median and 90th percentile,
w = 0.5 m/s) for ARG, FN and MBN on three configurations.
Costs: per-state CPU relative to six classes on the base grid, same host (production ratios for
NC1, NC11, NC9BC and NC6 x8)."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

CFG, LONG, OUT = sys.argv[1:4]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 6.5, "legend.frameon": False, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})
C = pd.read_csv(CFG)
PROD_COST = {"NC1@U1": 0.11, "NC11@U1": 1.82, "NC9BC@U1": 1.59, "NC6@U8": 87.9}
C["cost"] = [PROD_COST.get(c, v) for c, v in zip(C.cfg, C.cost_rel)]
FAM = {"NC1": "one class", "NC6": "6 classes (production groups)", "NC11": "6 groups, finer bounds",
       "NC9BC": "BC-resolved", "G1": "6 classes, salts and SOA regrouped", "G2": "hygroscopicity groups",
       "G3": "hygroscopicity groups + BC", "G4": "hygroscopicity + sea salt", "G5": "hygroscopicity + sea salt + BC"}
COL = {"one class": "#9e9e9e", "6 classes (production groups)": "#d9622b", "6 groups, finer bounds": "#e8a33d",
       "BC-resolved": "#7a5195", "6 classes, salts and SOA regrouped": "#b5523b", "hygroscopicity groups": "#3b7bb5",
       "hygroscopicity groups + BC": "#5fa8d3", "hygroscopicity + sea salt": "#1f7a5a",
       "hygroscopicity + sea salt + BC": "#7fbf9f"}
MK = {"U1": "o", "L2": "s", "L4": "^", "U8": "D", "U2": "v"}


def fam(cfg):
    s = cfg.split("@")[0]
    if s.startswith("G"):
        return FAM[s.split("_")[0]]
    return FAM[s]


def scatter(ax, sub, col, title, ylim):
    d = C[(C.subset == sub) & C[col].notna() & C.cost.notna()].copy()
    d = d[d.cfg.str.split("@").str[1].isin(MK)]
    for _, r in d.iterrows():
        f = fam(r.cfg); g = r.cfg.split("@")[1]
        ax.scatter(r.cost, r[col], marker=MK[g], s=26, color=COL[f], edgecolor="k", linewidth=0.3, zorder=3)
    # Pareto front (lower-left envelope)
    d = d.sort_values("cost"); best = np.inf; px, py = [], []
    for _, r in d.iterrows():
        if r[col] < best:
            best = r[col]; px.append(r.cost); py.append(r[col])
    ax.step(px, py, where="post", color="k", lw=0.7, ls="--", zorder=2)
    ax.set_xscale("log"); ax.set_xlim(0.08, 130); ax.set_ylim(*ylim)
    ax.set_xticks([0.1, 0.3, 1, 3, 10, 30, 100]); ax.set_xticklabels(["0.1", "0.3", "1", "3", "10", "30", "100"])
    ax.set_xlabel("Relative box-model cost")
    ax.set_title(title, fontsize=7.5, loc="left")
    return d


fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.0), gridspec_kw={"width_ratios": [1, 1, 1.05]})
scatter(axs[0], "strata", "S0.1_p90", "CCN(0.1%), all states, P90", (0, 40))
axs[0].set_ylabel("|simulation error| (%)")
scatter(axs[1], "strata_lowcloud", "S0.1_p90", "CCN(0.1%), low-cloud weighted, P90", (0, 20))
# (c) activation: end-to-end N_d, w = 0.5
ax = axs[2]
cfgs = ["NC6@U1", "NC6@L2", "G4_4@U1"]
cfgs = [c for c in cfgs if ((C.cfg == c) & (C.subset == "strata")).any()]
lab = {"arg": "ARG", "fn": "FN", "mbn": "FN, revised splitting"}
colk = {"arg": "#8c8c8c", "fn": "#7fbf9f", "mbn": "#1f7a5a"}
x0 = np.arange(len(cfgs)); wbar = 0.26
for k, (sch) in enumerate(["arg", "fn", "mbn"]):
    q = [C[(C.cfg == c) & (C.subset == "strata")][f"E_{sch}0.5_med"].iloc[0] for c in cfgs]
    p = [C[(C.cfg == c) & (C.subset == "strata")][f"E_{sch}0.5_p90"].iloc[0] for c in cfgs]
    ax.bar(x0 + (k - 1) * wbar, q, wbar * 0.9, color=colk[sch], label=lab[sch])
    ax.errorbar(x0 + (k - 1) * wbar, q, yerr=[np.zeros(len(q)), np.array(p) - np.array(q)], fmt="none", ecolor="k", elinewidth=0.6, capsize=1.5)
qa = [C[(C.cfg == c) & (C.subset == "strata")]["Sd0.5_med"].iloc[0] for c in cfgs]
ax.scatter(x0, qa, marker="_", s=120, color="k", zorder=4, label="aerosol only (parcel)")
NAME = {"NC6@U1": "6 classes\nbase grid", "NC6@L2": "6 classes\n2x 0.03-0.3 µm", "G4_4@U1": "4 classes\nbase grid"}
ax.set_xticks(x0); ax.set_xticklabels([NAME.get(c, c) for c in cfgs], fontsize=6.5)
ax.set_ylabel("|N$_d$ error|, w = 0.5 m s$^{-1}$ (%)"); ax.set_ylim(0, 100)
ax.set_title("N$_d$ end to end, median (whisker P90)", fontsize=7.5, loc="left")
ax.legend(loc="upper right", fontsize=6, bbox_to_anchor=(1.0, 1.0))
fams = sorted({fam(c) for c in C.cfg if c.split("@")[1] in MK}, key=list(COL).index)
h = [Line2D([], [], marker="o", ls="", color=COL[f], markeredgecolor="k", markeredgewidth=0.3, label=f) for f in fams]
h += [Line2D([], [], marker=MK[g], ls="", color="w", markeredgecolor="k", label={"U1": "base grid", "L2": "2x at 0.03-0.3 µm",
      "L4": "4x at 0.03-0.3 µm", "U8": "8x everywhere", "U2": "2x everywhere"}[g]) for g in ("U1", "L2", "L4", "U8")]
fig.legend(handles=h, loc="lower center", ncol=4, fontsize=6, bbox_to_anchor=(0.5, -0.01), columnspacing=0.8, handletextpad=0.2)
for a, t in zip(axs, "abc"):
    a.text(-0.16, 1.04, t, transform=a.transAxes, fontweight="bold", fontsize=9)
fig.tight_layout(rect=(0, 0.17, 1, 1), w_pad=1.2)
fig.savefig(f"{OUT}/fig4_cost_accuracy.png"); fig.savefig(f"{OUT}/fig4_cost_accuracy.pdf")
print("done")
