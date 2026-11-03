"""Process runs for candidate configurations: SCRAM box simulation of a scenario with a
configuration (class scheme @ size grid), end-state CCN, optics, N2O5 uptake, droplet number
(parcel / ARG / FN) and CPU time; time-mean |CCN error| against the particle-resolved reference.

  python tools/opt_run.py <scenario_id> config [config ...]
  env: OPT_SRC (runs/strata: PartMC reference), OPT_OUT (runs/opt/runs), OPT_SCRATCH
Writes OPT_OUT/<id>.json incrementally (resumable); the SCRAM run directory is removed."""
import copy, json, os, shutil, sys, tempfile, time
from pathlib import Path

import numpy as np

T_ = Path(__file__).resolve().parent; ROOT = T_.parent
sys.path.insert(0, str(T_)); sys.path.insert(0, str(ROOT))
from scrambox import scenarios as scn                          # noqa: E402
from scrambox import case as casemod                           # noqa: E402
from scrambox.post import run_all                              # noqa: E402
from scrambox.run import load_run, run_case                    # noqa: E402
import opt_configs as OC                                       # noqa: E402
import opt_screen as S                                         # noqa: E402  (fn_arg, parcel_on)

sid = sys.argv[1]
SRC = Path(os.environ.get("OPT_SRC", str(ROOT / "runs" / "strata"))) / sid
OUT = Path(os.environ.get("OPT_OUT", str(ROOT / "runs" / "opt" / "runs"))); OUT.mkdir(parents=True, exist_ok=True)
SCR = Path(os.environ.get("OPT_SCRATCH", str(Path(tempfile.gettempdir()) / "opt_runs"))) / sid
POST = {"ccn": {"enabled": True, "every": 2}, "optics": {"enabled": True, "every": 24},
        "hetchem": {"enabled": True, "every": 24, "gases": ["N2O5", "HO2"]}}


def base_case():
    c = copy.deepcopy(scn.get(sid)["case"])
    c["postprocess"] = POST
    if c["kernel"].get("isorropia", "v23") in ("none", None):
        c["kernel"]["isorropia"] = "v23"
    c["kernel"]["isorropia"] = c["kernel"].get("isorropia", "v23")
    return c


def main():
    tr = json.loads((SRC / "partmc_s1" / "truth.json").read_text())
    old = OUT / f"{sid}.json"                     # earlier per-scenario file (first batch)
    done = json.loads(old.read_text())["configs"] if old.exists() else {}
    (OUT / sid).mkdir(parents=True, exist_ok=True)
    for cfg in sys.argv[2:]:
        fo = OUT / sid / f"{cfg}.json"
        if fo.exists() or (cfg in done and "error" not in done[cfg]):
            continue
        res = {"id": sid, "configs": {}}
        mixing, ov = OC.mixing_and_overrides(cfg)
        c = casemod.set_path(base_case(), "mixing", mixing)
        for p, v in ov.items():
            c = casemod.set_path(c, p, v)
        d = SCR / cfg.replace("@", "_")
        shutil.rmtree(d, ignore_errors=True)
        t0 = time.time()
        try:
            run_case(c, d)
            run = load_run(d)
            summ = json.loads((d / "summary.json").read_text())
            post = json.loads((d / "post.json").read_text()) if (d / "post.json").exists() else run_all(run, POST, d)
            res_ccn = post["ccn"]["resolved"]
            ccn_t = {str(s["S_pct"]): [st["spectrum"][k]["ccn_cm3"] for st in res_ccn]
                     for k, s in enumerate(res_ccn[0]["spectrum"])}
            tpost = np.array(post["ccn"]["time_s"])
            m = {"ccn": {k: v[-1] for k, v in ccn_t.items()}, "cpu_s": summ.get("cpu_s"), "wall_s": None,
                 "nsize": int(run["grid"]["nsize"]), "nclass": int(run["grid"]["nfrac"]),
                 "returncode": summ.get("returncode")}
            # time-mean |relative error| of CCN against the reference at the reference output times
            tt = np.array(tr["time_s"]); tm = {}
            for s_, v in ccn_t.items():
                if s_ in tr["ccn_cm3"]:
                    sim = np.interp(tt, tpost, np.array(v)); ref = np.array(tr["ccn_cm3"][s_])
                    ok = (ref > 1) & (tt > 0)
                    tm[s_] = float(np.mean(np.abs(sim[ok] / ref[ok] - 1))) if ok.any() else None
            m["time_mean_abs"] = tm
            m["eabs"] = post["optics"]["resolved"][-1]["bc_abs_enhancement"]
            rr = post["hetchem"]["states"][-1]["rules"]["N2O5"]["resolved"]
            m["k_bt"] = (rr.get("bertram_thornton") or {}).get("k_s1"); m["k_cesm"] = rr["cesm"]["k_s1"]
            m["act"] = {str(w): {**S.parcel_on(d, w), **S.fn_arg(d, w)} for w in S.W}
            m["N_cm3"] = float(run["results"]["number_m3"][-1].sum() / 1e6)
        except Exception as e:
            m = {"error": repr(e)[:400]}
        m["wall_s"] = time.time() - t0
        res["configs"][cfg] = m
        tmp = fo.with_suffix(".tmp"); tmp.write_text(json.dumps(res)); os.replace(tmp, fo)
        shutil.rmtree(d, ignore_errors=True)
        print(sid, cfg, m.get("cpu_s"), round(m["wall_s"], 1), m.get("error", ""), flush=True)


if __name__ == "__main__":
    main()
