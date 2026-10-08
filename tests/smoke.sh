#!/bin/bash
# Smoke test: one case, three kernels (no-ISO, v1, v2.3).
cd /data/shupengz/Model/scram-boxapp2
python3 - <<'EOF'
import json, copy, time, numpy as np
from scrambox import case as cm
from scrambox.run import run_case
base = cm.load("tests/smoke_case.json")
for iso in ["none", "v1", "v23"]:
    c = copy.deepcopy(base); c["kernel"]["isorropia"] = iso
    if iso != "none":
        c["initial"]["gas_ugm3"].update({"NH3": 3.0, "HNO3": 2.0})
    try:
        out = run_case(c, f"runs/smoke_{iso}")
    except Exception as e:
        print(f"== {iso}: FAILED: {e}"); continue
    r, s, g = out["results"], out["summary"], out["grid"]
    sp = g["species"]["aerosol"]
    m = r["mass_ugm3"]; n = r["number_m3"]
    tot0, tot1 = m[0].sum(axis=(0, 1)), m[-1].sum(axis=(0, 1))
    print(f"== {iso}: substep={out['case']['time']['substep_s']} s, cpu={s['cpu_s']:.2f} s, wall={s['wall_s']:.1f} s, outputs={s['n_output']}, classes={g['nfrac']}")
    print("   N total (#/cm3): %.1f -> %.1f" % (n[0].sum()/1e6, n[-1].sum()/1e6))
    for i, name in enumerate(sp):
        if tot0[i] > 0 or tot1[i] > 0:
            print(f"   {name:7s} {tot0[i]:10.4f} -> {tot1[i]:10.4f} ug/m3")
    print("   SULF gas: %.4f -> %.4f" % (r["gas_ugm3"][0, g["species"]["gases"].index("SULF")], r["gas_ugm3"][-1, g["species"]["gases"].index("SULF")]))
    print("   conservation:", {k: f"{v:+.2e}" for k, v in s["conservation_rel"].items()})
    print("   notes:", s["notes"])
EOF
