#!/bin/bash
# Pure coagulation (no condensation, no gas) vs sub-step; internal mixture (NC1) and NC6.
cd /data/shupengz/Model/scram-boxapp2
python3 - <<'EOF'
import copy
from scrambox.run import run_case
base = {"case": {"name": "coag only"}, "kernel": {"isorropia": "none"},
        "time": {"duration_s": 3600, "dt_s": 600, "output_every_s": 600},
        "processes": {"coagulation": True, "condensation": False},
        "initial": {"particles": {"mode": "lognormal", "modes": [
            {"name": "acc", "number_cm3": 3000, "dg_um": 0.10, "sigma_g": 1.6, "mass_fraction": {"PSO4": 1.0}},
            {"name": "ait", "number_cm3": 2000, "dg_um": 0.05, "sigma_g": 1.7, "mass_fraction": {"PSO4": 1.0}}]}}}
print("mix   substep   loss%(1h)   per-call loss (#/cm3)   dry mass change")
for mix in ["NC1", "NC6"]:
    for sub in [600, 300, 60, 10, 1]:
        c = copy.deepcopy(base); c["mixing"] = {"preset": mix}; c["time"]["substep_s"] = sub
        out = run_case(c, f"runs/coag_{mix}_{sub}")
        r = out["results"]; n = r["number_m3"]
        iw = out["grid"]["species"]["aerosol"].index("PH2O")
        m = r["mass_ugm3"].sum(axis=(1, 2)); m = m.sum(axis=1) - m[:, iw]   # dry: water is diagnosed
        loss = n[0].sum() - n[-1].sum()
        ncall = 3600 / sub
        print(f"{mix:4s}  {sub:6d}   {100*loss/n[0].sum():7.2f}     {loss/1e6/ncall:12.4f}          {(m[-1]-m[0])/m[0]:+.1e}")
EOF
