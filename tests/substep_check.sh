#!/bin/bash
# Is the number loss driven by the sub-step or by the ISORROPIA path? Coagulation on/off too.
cd ${BOX_REPO:?}
python3 - <<'EOF'
import copy
from scrambox import case as cm
from scrambox.run import run_case
base = cm.load("tests/smoke_case.json")
rows = []
for iso in ["none", "v23"]:
    for sub in [600, 60, 10]:
        for coag in [True, False]:
            c = copy.deepcopy(base); c["kernel"]["isorropia"] = iso
            c["time"]["substep_s"] = sub; c["processes"] = {"coagulation": coag}
            if iso != "none":
                c["initial"]["gas_ugm3"].update({"NH3": 3.0, "HNO3": 2.0})
            out = run_case(c, f"runs/sub_{iso}_{sub}_{int(coag)}")
            n = out["results"]["number_m3"]
            rows.append((iso, sub, coag, n[0].sum()/1e6, n[-1].sum()/1e6, out["summary"]["cpu_s"]))
print("isorropia substep coag   N0       N(1h)    loss%   cpu")
for r in rows:
    print(f"{r[0]:9s} {r[1]:5d}   {str(r[2]):5s} {r[3]:8.1f} {r[4]:8.1f} {100*(1-r[4]/r[3]):6.2f}  {r[5]:.2f}")
EOF
