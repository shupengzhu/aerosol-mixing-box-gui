
import os, glob, json, shutil, tarfile, subprocess, sys
import numpy as np, pandas as pd
def merge_and_analyse(t118, t86, root, ss="strata_summary_126.json", st="strata15.json"):
    shutil.rmtree(root, ignore_errors=True); os.makedirs(f"{root}/h118"); os.makedirs(f"{root}/h86")
    with tarfile.open(t118) as t: t.extractall(f"{root}/h118")
    with tarfile.open(t86) as t: t.extractall(f"{root}/h86")
    M = f"{root}/runs"; shutil.copytree(f"{root}/h118/runs", M)
    for f in glob.glob(f"{root}/h86/out/*/*.json"):
        sid = os.path.basename(os.path.dirname(f)); os.makedirs(f"{M}/{sid}", exist_ok=True)
        dst = f"{M}/{sid}/{os.path.basename(f)}"
        if not os.path.exists(dst): shutil.copy(f, dst)
    cp = []
    def items(rt, h):
        for f in glob.glob(f"{rt}/*.json") + glob.glob(f"{rt}/*/*.json"):
            d = json.load(open(f))
            for c, m in d["configs"].items(): cp.append((d["id"], c, h, m.get("cpu_s")))
    items(f"{root}/h118/runs", "118"); items(f"{root}/h86/out", "86"); items(f"{root}/h86/out_norm", "86n")
    CP = pd.DataFrame(cp, columns=["id", "cfg", "host", "cpu"]); CP["hg"] = CP.host.str.replace("86n", "86")
    base = CP[CP.cfg == "NC6@U1"].groupby(["id", "hg"]).cpu.median()
    CP["cost_rel"] = [c / base.get((i, h), np.nan) for i, h, c in zip(CP.id, CP.hg, CP.cpu)]
    CP[CP.host != "86n"].groupby(["id", "cfg"]).cost_rel.median().reset_index().to_csv(f"{root}/cost_rel.csv", index=False)
    r = subprocess.run([sys.executable, "grl/opt/opt_analysis.py", f"{root}/h118/screen", M, ss, "regimes_summary.json", st, f"{root}/an", f"{root}/cost_rel.csv"], capture_output=True, text=True)
    return r.returncode, r.stderr[-800:]
