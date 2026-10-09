"""WP2 benchmark: SCRAM box (NC1 / NC6 / NC11) against the PartMC particle-resolved
reference driven by the same scenario.

  python3 tools/partmc_experiment.py <scenario_id> [--n-part 20000] [--seeds 1 2 3] [--out runs/wp2]

For every scenario:
  scram_<NC>            SCRAM box runs, CESM post-processing (optics, CCN, het. chemistry)
  partmc_s<seed>        PartMC reference runs (pmc_scram_thermo), particle-resolved truth.json
  partmc_s<seed>_<NC>   the PartMC particles projected on the same NC grid, same post-processing
The differences separate (a) the error of SCRAM's sectional/class dynamics
(scram_NC vs partmc_NC: same representation, different evolution) from (b) the
representation (discretisation) error (partmc_NC vs particle-resolved truth).
Writes summary.json and PNG figures in the output directory.
"""
from __future__ import annotations

import argparse
import copy
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scrambox import scenarios as scn                      # noqa: E402
from scrambox.partmc import chi_sections, project, run_partmc, truth   # noqa: E402
from scrambox.post import run_all                          # noqa: E402
from scrambox.run import load_run, run_case                # noqa: E402

PRESETS = ["NC1", "NC6", "NC11"]
POST = {"optics": {"enabled": True, "every": 2}, "ccn": {"enabled": True, "every": 2},
        "hetchem": {"enabled": True, "every": 2, "gases": ["N2O5", "HO2"]}}


def series(run):
    r, g = run["results"], run["grid"]
    aer = g["species"]["aerosol"]
    m = r["mass_ugm3"].sum(axis=(1, 2))
    iw = aer.index("PH2O")
    out = {"time_h": (r["time_s"] / 3600).tolist(),
           "N_cm3": (r["number_m3"].sum(axis=(1, 2)) / 1e6).tolist(),
           "PM_dry": (m.sum(axis=1) - m[:, iw]).tolist(),
           "species": {s: m[:, i].tolist() for i, s in enumerate(aer) if m[:, i].max() > 0},
           "chi_sections": [chi_sections(r["mass_ugm3"][k], r["number_m3"][k], aer) for k in range(len(r["time_s"]))]}
    post = run.get("post") or {}
    if "ccn" in post:
        c = post["ccn"]
        out["post_time_h"] = [t / 3600 for t in c["time_s"]]
        out["ccn"] = {str(s["S_pct"]): [st["spectrum"][k]["ccn_cm3"] for st in c["resolved"]]
                      for k, s in enumerate(c["resolved"][0]["spectrum"])}
        out["nd"] = {str(a["w_m_s"]): [st["activation"][k]["nd_cm3"] for st in c["resolved"]]
                     for k, a in enumerate(c["resolved"][0]["activation"])}
    if "optics" in post:
        o = post["optics"]
        out["ext550"] = [x["ext550_Mm"] for x in o["resolved"]]
        out["ssa550"] = [x["ssa550"] for x in o["resolved"]]
        out["abs550"] = [x["abs550_Mm"] for x in o["resolved"]]
        out["bc_abs_enh"] = [x["bc_abs_enhancement"] for x in o["resolved"]]
    if "hetchem" in post:
        h = post["hetchem"]
        out["k_n2o5_bt"] = [st["rules"]["N2O5"]["resolved"].get("bertram_thornton", {}).get("k_s1") for st in h["states"]]
        out["k_n2o5_cesm"] = [st["rules"]["N2O5"]["resolved"]["cesm"]["k_s1"] for st in h["states"]]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("scenario")
    ap.add_argument("--n-part", type=int, default=20000)
    ap.add_argument("--seeds", type=int, nargs="+", default=[1, 2, 3])
    ap.add_argument("--out", default="runs/wp2")
    ap.add_argument("--hours", type=float, default=None)
    ap.add_argument("--report-only", action="store_true", help="rebuild report.md from summary.json")
    ap.add_argument("--post-only", action="store_true",
                    help="reuse finished SCRAM and PartMC runs; redo projections, truth and post-processing")
    a = ap.parse_args()
    out = Path(a.out) / a.scenario
    out.mkdir(parents=True, exist_ok=True)
    if a.report_only:
        report(out, json.loads((out / "summary.json").read_text()))
        return
    base = copy.deepcopy(scn.get(a.scenario)["case"])
    if a.hours:
        base["time"]["duration_s"] = a.hours * 3600
    base["postprocess"] = POST
    base["kernel"]["isorropia"] = base["kernel"].get("isorropia", "v23")
    if base["kernel"]["isorropia"] == "none":
        base["kernel"]["isorropia"] = "v23"

    def scram(nc):
        if a.post_only:
            return nc
        c = copy.deepcopy(base); c["mixing"] = {"preset": nc}
        run_case(c, out / f"scram_{nc}")
        return nc

    def partmc(seed):
        rd = out / f"partmc_s{seed}"
        if not a.post_only:
            run_partmc(base, rd, n_part=a.n_part, seed=seed)
        else:   # PartMC output exists; redo only the steps after the run that are missing
            if not (rd / "summary.json").exists():
                info = project(rd)
                info["summary"].update(n_part=a.n_part, model="partmc")
                (rd / "summary.json").write_text(json.dumps(info["summary"], indent=1))
            if not (rd / "truth.json").exists():
                truth(rd)
        for nc in PRESETS:
            d = out / f"partmc_s{seed}_{nc}"
            project(out / f"partmc_s{seed}", d, {"preset": nc})
            run_all(load_run(d), POST, d)
        return seed

    with ThreadPoolExecutor(max_workers=len(PRESETS) + len(a.seeds)) as ex:
        jobs = [ex.submit(scram, nc) for nc in PRESETS] + [ex.submit(partmc, s) for s in a.seeds]
        for j in jobs:
            print("done", j.result(), flush=True)

    summary = {"scenario": a.scenario, "n_part": a.n_part, "seeds": a.seeds, "runs": {}, "truth": {}}
    for nc in PRESETS:
        summary["runs"][f"scram_{nc}"] = series(load_run(out / f"scram_{nc}"))
        for s in a.seeds:
            summary["runs"][f"partmc_s{s}_{nc}"] = series(load_run(out / f"partmc_s{s}_{nc}"))
    for s in a.seeds:
        summary["truth"][f"s{s}"] = json.loads((out / f"partmc_s{s}" / "truth.json").read_text())
    (out / "summary.json").write_text(json.dumps(summary))
    report(out, summary)


def _mean_sd(summary, key_fn, seeds):
    arr = np.array([key_fn(s) for s in seeds], dtype=float)
    return arr.mean(axis=0), arr.std(axis=0)


def report(out, S):
    seeds = S["seeds"]
    R, T = S["runs"], S["truth"]
    lines = [f"# WP2 benchmark: {S['scenario']}", "",
             f"PartMC reference: {S['n_part']} computational particles, seeds {seeds}; "
             "mean ± sd over seeds.", ""]
    th = np.array(T[f"s{seeds[0]}"]["time_s"]) / 3600
    end = -1
    lines.append("| quantity (end of run) | " + " | ".join(f"SCRAM {nc}" for nc in PRESETS) + " | "
                 + " | ".join(f"PartMC→{nc}" for nc in PRESETS) + " | PartMC truth |")
    lines.append("|---|" + "---|" * (2 * len(PRESETS) + 1))

    def f(x):
        return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else fmt.format(x)

    fmt = "{:.3g}"

    def ms(vals):
        v = [x for x in vals if x is not None]
        return (f(float(np.mean(v))) + " ± " + f(float(np.std(v)))) if v else "–"

    def row(label, key, sub=None, truth=None):
        if key not in R[f"scram_{PRESETS[0]}"]:
            return
        if sub is not None and any(sub not in R[k][key] for k in R):
            return
        cells = []
        for nc in PRESETS:
            v = R[f"scram_{nc}"][key] if sub is None else R[f"scram_{nc}"][key][sub]
            cells.append(f(v[end]))
        for nc in PRESETS:
            cells.append(ms([(R[f"partmc_s{s}_{nc}"][key] if sub is None else R[f"partmc_s{s}_{nc}"][key][sub])[end]
                             for s in seeds]))
        tcell = ms([truth(T[f"s{s}"]) for s in seeds]) if truth else "–"
        lines.append(f"| {label} | " + " | ".join(cells) + f" | {tcell} |")

    row("N (cm⁻³)", "N_cm3", truth=lambda t: t["n_total_cm3"][-1])
    row("PM dry (µg m⁻³)", "PM_dry")
    for sp in ["PBC", "PSO4", "PNH4", "PNO3", "PHCL", "PSOA3P", "PMD"]:
        if sp in R[f"scram_{PRESETS[0]}"]["species"]:
            row(f"{sp} (µg m⁻³)", "species", sp)
    row("χ mixing-state index", "chi_sections", truth=lambda t: t["chi"][-1])
    for s_ in ["0.1", "0.2", "0.5"]:
        row(f"CCN S={s_}% (cm⁻³)", "ccn", s_, truth=lambda t, s_=s_: t["ccn_cm3"][s_][-1])
    for w in ["0.2", "1.0"]:
        row(f"Nd w={w} m/s (cm⁻³)", "nd", w)
    row("ext 550 nm (Mm⁻¹)", "ext550")
    row("abs 550 nm (Mm⁻¹)", "abs550")
    row("BC E_abs", "bc_abs_enh")
    row("k N2O5 BT09 (s⁻¹)", "k_n2o5_bt")
    lines += ["", "Columns: SCRAM box evolution in each representation; PartMC particles projected on the same "
              "representation and processed with the same CESM post-processing; PartMC particle-resolved values "
              "(CCN from per-particle κ-Köhler, χ from per-particle composition).", ""]
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        figs = [("N_cm3", None, "N (cm⁻³)", lambda t: t["n_total_cm3"]),
                ("chi_sections", None, "mixing-state index χ", lambda t: t["chi"]),
                ("ccn", "0.2", "CCN S=0.2 % (cm⁻³)", lambda t: t["ccn_cm3"]["0.2"]),
                ("ccn", "0.1", "CCN S=0.1 % (cm⁻³)", lambda t: t["ccn_cm3"]["0.1"]),
                ("nd", "0.2", "Nd w=0.2 m/s (cm⁻³)", None),
                ("bc_abs_enh", None, "BC absorption enhancement", None),
                ("species", "PNO3", "nitrate (µg m⁻³)", None),
                ("k_n2o5_bt", None, "k(N2O5), BT09 (s⁻¹)", None)]
        fig, axs = plt.subplots(2, 4, figsize=(18, 8))
        colors = {"NC1": "tab:blue", "NC6": "tab:orange", "NC11": "tab:green"}
        for ax, (key, sub, title, tf) in zip(axs.ravel(), figs):
            for nc in PRESETS:
                r = R[f"scram_{nc}"]
                if key not in r or (sub and sub not in r[key]):
                    continue
                x = r["post_time_h"] if key in ("ccn", "nd", "bc_abs_enh", "k_n2o5_bt", "ext550") else r["time_h"]
                y = r[key][sub] if sub else r[key]
                ax.plot(x, y, color=colors[nc], lw=2, label=f"SCRAM {nc}")
                ys = [np.array((R[f"partmc_s{s}_{nc}"][key][sub] if sub else R[f"partmc_s{s}_{nc}"][key]), dtype=float) for s in seeds]
                xs = R[f"partmc_s{seeds[0]}_{nc}"]["post_time_h"] if key in ("ccn", "nd", "bc_abs_enh", "k_n2o5_bt", "ext550") else R[f"partmc_s{seeds[0]}_{nc}"]["time_h"]
                ax.plot(xs, np.mean(ys, axis=0), color=colors[nc], ls="--", lw=1.5, label=f"PartMC→{nc}")
            if tf:
                ys = [np.array(tf(T[f"s{s}"]), dtype=float) for s in seeds]
                ax.plot(th, np.mean(ys, axis=0), color="k", lw=2.2, label="PartMC truth")
            ax.set_title(title); ax.set_xlabel("time (h)")
        axs[0, 0].legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(out / "benchmark.png", dpi=110)
        lines.append("![benchmark](benchmark.png)")
    except ImportError:
        lines.append("(matplotlib not available: no figure)")
    (out / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
