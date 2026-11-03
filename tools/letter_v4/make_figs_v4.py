"""Letter v4, Figures 1 and 3.
python make_figs_v4.py <opt_long.csv> <ensemble_table.csv> <mam4_ndarg.json> <synth732.json> <out_dir>
Fig 1  the error chain by region type: (a) CCN(0.1%): MAM4 modes, six-class representation,
       six-class simulation, recommended configuration; (b) N_d at w = 0.5 m/s: MAM4 modes with the
       parcel model and with ARG, six-class simulation with the parcel model, ARG, FN and FN with the
       revised population splitting.
Fig 3  activation: (a) activation factor by scheme and FN branch on the simulated states, three
       updrafts; (b) FN and revised splitting on 732 synthetic populations by branch."""
import json, math, os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.ticker
import matplotlib.pyplot as plt

LONG, ENS, MND, SYN, OUT = sys.argv[1:6]; os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "font.family": "DejaVu Sans", "axes.linewidth": 0.6, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "legend.fontsize": 7, "legend.frameon": False, "savefig.dpi": 300,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = {"mam4": "#9e9e9e", "mam4arg": "#5e5e5e", "rep": "#3b7bb5", "sim": "#d9622b", "best": "#f3b48a",
       "arg": "#8c8c8c", "fn": "#7fbf9f", "mbn": "#1f7a5a"}
L = pd.read_csv(LONG); E = pd.read_csv(ENS); M = json.load(open(MND))
REC = "NC6@L2"
base = L[L.cfg == "NC6@U1"].set_index("id"); best = L[L.cfg == REC].set_index("id")
T = base[["layer", "type"]].copy()
for c in ("R0.1", "S0.1", "Sd0.5", "E_arg0.5", "E_fn0.5", "E_mbn0.5", "nd_t0.5"):
    T[c] = base[c]
T["B0.1"] = best["S0.1"]
e = E[~E.internal.astype(bool)].set_index("id")
T["M0.1"] = e["M0.1"]; T["Mp0.5"] = e["A_mam0.5"]
T["Ma0.5"] = [((M.get(i) or {}).get("nd_arg") or {}).get("0.5", np.nan) / t if t and t > 1 else np.nan for i, t in zip(T.index, T["nd_t0.5"])]
order = ["polluted land", "biomass burning", "dust", "ocean", "remote", "global background"]
labels = {"polluted land": "polluted land\n(n=%d)", "biomass burning": "biomass\nburning (n=%d)", "dust": "dust (n=%d)",
          "ocean": "ocean (n=%d)", "remote": "remote (n=%d)", "global background": "global\nbackground (n=%d)"}
pct = lambda x: 100 * (np.asarray(x, float) - 1)
rng = np.random.default_rng(3)


def strip(ax, D, cols, ylab=True, xlim=(-90, 250)):
    for ti, typ in enumerate(order):
        sel = D[D.type == typ]
        for j, (k, col, lab) in enumerate(cols):
            v = pct(sel[k]); v = v[np.isfinite(v)]
            if len(v) == 0:
                continue
            y = ti * (len(cols) + 1) + j
            ax.scatter(v, y + rng.uniform(-0.28, 0.28, len(v)), s=4, color=col, alpha=0.6, lw=0, label=lab if ti == 0 else None, zorder=2)
            ax.plot([np.median(v)] * 2, [y - 0.45, y + 0.45], color="k", lw=1.2, zorder=3)
    ax.axvline(0, color="k", lw=0.6); ax.axvspan(-5, 5, color="0.92", zorder=0, lw=0)
    ax.set_yticks([ti * (len(cols) + 1) + (len(cols) - 1) / 2 for ti in range(len(order))])
    ax.set_yticklabels([labels[t] % len(D[D.type == t]) for t in order] if ylab else [""] * len(order))
    ax.invert_yaxis(); ax.set_xscale("symlog", linthresh=10, linscale=1.0); ax.set_xlim(*xlim)
    ax.set_xticks([-50, -10, 0, 10, 50, 200]); ax.set_xticklabels(["-50", "-10", "0", "10", "50", "200"]); ax.tick_params(axis="y", length=0)


fig, axs = plt.subplots(1, 2, figsize=(7.2, 4.4))
ca = [("M0.1", COL["mam4"], "MAM4 modes"), ("R0.1", COL["rep"], "six classes, particles projected"),
      ("S0.1", COL["sim"], "six classes, simulated"), ("B0.1", COL["best"], "six classes, 2x sections 0.03-0.3 µm")]
strip(axs[0], T, ca); axs[0].set_title("CCN at 0.1% supersaturation", fontsize=8, loc="left")
cb = [("Ma0.5", COL["mam4arg"], "MAM4 modes + ARG"), ("Sd0.5", COL["sim"], "six classes, simulated + parcel"),
      ("E_fn0.5", COL["fn"], "six classes, simulated + FN"), ("E_mbn0.5", COL["mbn"], "six classes, simulated + revised FN")]
strip(axs[1], T, cb, ylab=False, xlim=(-95, 250)); axs[1].set_title("Droplet number, w = 0.5 m s$^{-1}$", fontsize=8, loc="left")
for ax in axs:
    ax.set_xlabel("Error relative to particle-resolved reference (%)", fontsize=7)
fig.tight_layout(rect=(0, 0.22, 1, 0.97), w_pad=1.0)
for ax, t_ in zip(axs, "ab"):
    h, l = ax.get_legend_handles_labels(); bb = ax.get_position()
    fig.legend(h, l, loc="upper left", bbox_to_anchor=(bb.x0, 0.215), ncol=1, fontsize=6.3,
               markerscale=2.5, handletextpad=0.1, columnspacing=0.5)
    fig.text(bb.x0 - (0.13 if t_ == "a" else 0.02), bb.y1 + 0.02, t_, fontweight="bold", fontsize=9)
fig.savefig(f"{OUT}/fig1_error_chain.png"); fig.savefig(f"{OUT}/fig1_error_chain.pdf"); plt.close(fig)

# ---------------------------------------------------------------- Figure 3
S = L[(L.cfg == "NC6@U1") & (L.layer == "strata")]
fig, axs = plt.subplots(1, 2, figsize=(7.2, 2.8), gridspec_kw={"width_ratios": [1.5, 1]})
ax = axs[0]
groups = [("A_arg", None, COL["arg"], "ARG"), ("A_fn", True, COL["fn"], "FN, analytic branch"), ("A_fn", False, "#b9dfca", "FN, empirical branch"),
          ("A_mbn", True, COL["mbn"], "revised, analytic"), ("A_mbn", False, "#5aa37f", "revised, empirical")]
for wi, w in enumerate(("0.2", "0.5", "1.0")):
    D = S[f"Dpos{w}"].astype(str).str.lower() == "true"
    for gi, (k, br, col, lab) in enumerate(groups):
        v = S[f"{k}{w}"] if br is None else S[f"{k}{w}"][D if br else ~D]
        v = v[np.isfinite(v) & (v > 0)]
        x = wi * (len(groups) + 1.5) + gi
        ax.scatter(x + rng.uniform(-0.3, 0.3, len(v)), v, s=4, color=col, alpha=0.6, lw=0, label=lab if wi == 0 else None)
        ax.plot([x - 0.4, x + 0.4], [np.median(v)] * 2, color="k", lw=1.2)
ax.axhspan(0.9, 1.1, color="0.92", zorder=0, lw=0); ax.axhline(1, color="k", lw=0.6)
ax.set_yscale("log"); ax.set_ylim(0.1, 2.5); ax.set_yticks([0.125, 0.25, 0.5, 1, 2]); ax.set_yticklabels(["0.125", "0.25", "0.5", "1", "2"])
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks([wi * (len(groups) + 1.5) + 2 for wi in range(3)]); ax.set_xticklabels(["0.2", "0.5", "1"])
ax.set_xlabel("Updraft w (m s$^{-1}$)"); ax.set_ylabel("Parameterized / parcel N$_d$")
ax.set_title("Activation factor on the simulated states (six classes)", fontsize=7.5, loc="left")
ax.legend(loc="lower right", ncol=2, fontsize=6, markerscale=2.5, handletextpad=0.1)
ax = axs[1]
R = json.load(open(SYN))
K = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "e3_constants.json")))
al, G, A = K["alpha"], K["G"], K["A"]
rr = [r for r in R if r["truth"] > 1]
Dp = np.array([r["smax_FN"] ** 4 - 16 * A ** 2 * al * r["w"] / (9 * G) >= 0 for r in rr])
fn = np.array([r["FN"] / r["truth"] for r in rr]); mb = np.array([r["MBN"] / r["truth"] for r in rr]); ct = np.array([r["count_at_true_smax"] / r["truth"] for r in rr])
bins = np.geomspace(0.25, 1.6, 26)
ax.hist(fn[~Dp], bins, color="#b9dfca", label=f"FN, empirical (median {np.median(fn[~Dp]):.2f})")
ax.hist(mb[~Dp], bins, histtype="step", color=COL["mbn"], lw=1.3, label=f"revised (median {np.median(mb[~Dp]):.2f})")
ax.hist(ct[~Dp], bins, histtype="step", color="k", lw=0.8, ls="--", label=f"counting at exact S$_{{max}}$ ({np.median(ct[~Dp]):.2f})")
ax.set_xscale("log"); ax.set_xticks([0.25, 0.5, 1]); ax.set_xticklabels(["0.25", "0.5", "1"]); ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.axvline(1, color="k", lw=0.6); ax.set_xlabel("Parameterized / parcel N$_d$"); ax.set_ylabel("Synthetic cases")
ax.set_title(f"Empirical-branch cases (n={int((~Dp).sum())} of {len(rr)})", fontsize=7.5, loc="left")
ax.legend(loc="upper left", fontsize=6)
for a, t in zip(axs, "ab"):
    a.text(-0.12, 1.04, t, transform=a.transAxes, fontweight="bold", fontsize=9)
fig.tight_layout(); fig.savefig(f"{OUT}/fig3_activation.png"); fig.savefig(f"{OUT}/fig3_activation.pdf"); plt.close(fig)
print("done")
