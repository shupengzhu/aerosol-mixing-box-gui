#!/bin/bash
# Generate the scenario library, describe every scenario, then run each one (short
# version unless FULL=1) with its own kernel; the mixing-state scenarios with v1 and v23.
cd /data/shupengz/Model/scram-boxapp2
python3 tools/make_scenarios.py > /dev/null || exit 1
python3 - <<'EOF'
import copy, json, os, time
from scrambox import scenarios as sc
from scrambox.run import run_case
full = os.environ.get("FULL") == "1"
print(f"{'scenario':26s} {'kern':4s} {'N0 cm-3':>9s} {'PM10_0':>7s} {'N_end':>9s} {'PM10_end':>8s} {'NO3_end':>7s} {'cpu s':>7s}")
for s in sc.index():
    full_s = sc.get(s["id"])
    d = sc.describe(full_s["case"])
    kerns = ["v1", "v23"] if s["category"] == "mixing_state" else [full_s["case"]["kernel"]["isorropia"]]
    for k in kerns:
        c = copy.deepcopy(full_s["case"]); c["kernel"]["isorropia"] = k
        if not full:
            c["time"]["duration_s"] = min(c["time"]["duration_s"], 3600.0)
            c["time"]["output_every_s"] = 600.0
        out = run_case(c, f"runs/scn_{s['id']}_{k}")
        r, g = out["results"], out["grid"]
        aer = g["species"]["aerosol"]; iw = aer.index("PH2O")
        m = r["mass_ugm3"]; dry = m.sum(axis=(1, 2)); dry = dry.sum(axis=1) - dry[:, iw]
        no3 = m[-1, :, :, aer.index("PNO3")].sum() if "PNO3" in aer else 0.0
        print(f"{s['id']:26s} {k:4s} {r['number_m3'][0].sum()/1e6:9.1f} {dry[0]:7.2f} "
              f"{r['number_m3'][-1].sum()/1e6:9.1f} {dry[-1]:8.2f} {no3:7.3f} {out['summary']['cpu_s']:7.1f}")
    if d["notes"]:
        print("   notes:", d["notes"])
EOF
