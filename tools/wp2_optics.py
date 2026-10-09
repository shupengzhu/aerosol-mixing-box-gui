"""Optics part of the WP2 benchmark (run after tools/partmc_experiment.py):
particle-resolved optics of the PartMC runs against sectional optics of the SCRAM runs and
of the PartMC particles projected on NC1 / NC6 / NC11, all with the CESM reduced-shell
optics and the dry radius from each population's own N and M (prognostic moments).

  python3 tools/wp2_optics.py runs/wp2/<scenario> [--every 4]
Appends an optics table to report.md and writes optics.json.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scrambox.partmc import truth_optics          # noqa: E402
from scrambox.post import optics                  # noqa: E402
from scrambox.post.common import populations      # noqa: E402
from scrambox.run import load_run                 # noqa: E402

PRESETS = ["NC1", "NC6", "NC11"]
OPTS = {"mixing": 2, "perclass": True, "prognostic_moments": True, "dryshell_selector": 0}


def sectional(run_dir, times):
    r = load_run(run_dir)
    ts = r["results"]["time_s"]
    out = []
    for t in times:
        k = int(np.argmin(abs(ts - t)))
        res = optics.compute(populations(r, k), OPTS)
        out.append({key: res[key] for key in ("ext550_Mm", "abs550_Mm", "ssa550", "bc_abs_enhancement")})
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp_dir")
    ap.add_argument("--every", type=int, default=4)
    a = ap.parse_args()
    D = Path(a.exp_dir)
    seeds = sorted(int(p.name.split("_s")[1]) for p in D.glob("partmc_s*") if p.name.count("_") == 1)
    truth = {s: truth_optics(D / f"partmc_s{s}", every=a.every, opts=OPTS) for s in seeds}
    times = truth[seeds[0]]["time_s"]
    res = {"time_s": times, "truth": truth, "scram": {}, "partmc": {}}
    for nc in PRESETS:
        res["scram"][nc] = sectional(D / f"scram_{nc}", times)
        res["partmc"][nc] = {s: sectional(D / f"partmc_s{s}_{nc}", times) for s in seeds}
    (D / "optics.json").write_text(json.dumps(res))

    def ms(vals):
        v = np.array([x for x in vals if x is not None], dtype=float)
        return f"{v.mean():.3g} ± {v.std():.2g}" if v.size else "–"

    lines = ["", "## Optics (550 nm, CESM reduced-shell LUT, dry radius from N and M)", "",
             "| time (h) | quantity | " + " | ".join(f"SCRAM {nc}" for nc in PRESETS) + " | "
             + " | ".join(f"PartMC→{nc}" for nc in PRESETS) + " | PartMC particle-resolved |",
             "|---|---|" + "---|" * (2 * len(PRESETS) + 1)]
    for k in (len(times) // 2, len(times) - 1):
        for q, lab in (("ext550_Mm", "ext (Mm⁻¹)"), ("abs550_Mm", "abs (Mm⁻¹)"), ("bc_abs_enhancement", "BC E_abs")):
            cells = [f"{res['scram'][nc][k][q]:.3g}" if res['scram'][nc][k][q] is not None else "–" for nc in PRESETS]
            cells += [ms([res["partmc"][nc][s][k][q] for s in seeds]) for nc in PRESETS]
            cells.append(ms([truth[s][q][k] for s in seeds]))
            lines.append(f"| {times[k] / 3600:.0f} | {lab} | " + " | ".join(cells) + " |")
    with open(D / "report.md", "a") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
