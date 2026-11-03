"""Representation screening: PartMC end-state particles projected on every configuration
(class scheme @ size grid), CCN, parcel / ARG / FN droplet number, and (U1 grids only) optics
and N2O5 uptake.  No SCRAM run is needed.

  python tools/opt_screen.py <scenario_id> [config ...]
  env: OPT_SRC (default runs/strata), OPT_OUT (runs/opt/screen), OPT_SCRATCH (system temp dir)
Writes OPT_OUT/<id>.json incrementally (resumable)."""
import collections, json, os, shutil, sys, tempfile, time
from pathlib import Path

T_ = Path(__file__).resolve().parent; ROOT = T_.parent
sys.path.insert(0, str(T_)); sys.path.insert(0, str(ROOT))
import scrambox.partmc as PM                                  # noqa: E402
from scrambox.post import run_all                             # noqa: E402
from scrambox.post import ccn as C                            # noqa: E402
from scrambox.post.common import populations                  # noqa: E402
from scrambox.run import load_run                             # noqa: E402
from parcel_activation import parcel_nd                       # noqa: E402
from wp2_activation_parcel import scram_pops                  # noqa: E402
import activation_validation as av                            # noqa: E402
import opt_configs as OC                                      # noqa: E402

sid = sys.argv[1]
SRC = Path(os.environ.get("OPT_SRC", str(ROOT / "runs" / "strata"))) / sid
OUT = Path(os.environ.get("OPT_OUT", str(ROOT / "runs" / "opt" / "screen"))); OUT.mkdir(parents=True, exist_ok=True)
SCR = Path(os.environ.get("OPT_SCRATCH", str(Path(tempfile.gettempdir()) / "opt_screen"))) / sid
cfgs = sys.argv[2:] or OC.all_configs()
W = (0.2, 0.5, 1.0)
POST_C = {"ccn": {"enabled": True, "every": 1}}
POST_F = {"ccn": {"enabled": True, "every": 1}, "optics": {"enabled": True, "every": 1},
          "hetchem": {"enabled": True, "every": 1, "gases": ["N2O5", "HO2"]}}

_orig, _cache = PM.read_particles, {}


def _last(d):                       # project() only needs the end state: read the files once
    k = str(d)
    if k not in _cache:
        _cache[k] = collections.deque(_orig(d), maxlen=1)[0]
    yield _cache[k]


PM.read_particles = _last


def fn_arg(run_dir, w):
    pop = populations(load_run(run_dir), -1); md = C.modes(pop)
    av.T, av.P = pop.T, pop.P; zeta, pref, A = av.arg_consts(w)
    sl = [{"N": m["n_m3"], "sc": m["sc"], "dlo": m["dlo"], "dhi": m["dhi"], "tilt": m["tilt"], "kappa": m["kappa"],
           "dmid": m["dmid"]} for m in md]
    s_arg = av.maxsat([s["N"] for s in sl], [s["sc"] for s in sl], 0.5, 1.0, zeta, pref)
    Nf, scf = av.fn_spectrum(sl); st = []; s_fn = av.smax_fn(Nf, scf, w, pop.T, pop.P, stats=st)
    s_mbn = av.smax_mbn(Nf, scf, w, pop.T, pop.P)
    nd = lambda s: sum(C._activated(m, s) for m in md) / 1e6
    return {"arg": nd(s_arg), "arg_smax": s_arg * 100, "fn": nd(s_fn), "fn_smax": s_fn * 100, "fn_D_pos": bool(st[0]),
            "mbn": nd(s_mbn), "mbn_smax": s_mbn * 100}


def parcel_on(run_dir, w):
    pop, nn, rd, kk = scram_pops(run_dir); sm, nd = parcel_nd(pop.T, pop.P, w, nn, rd, kk)
    return {"parcel": nd / 1e6, "parcel_smax": sm * 100}


def metrics(d, full):
    run = load_run(d)
    post = run_all(run, POST_F if full else POST_C, None)
    out = {"ccn": {str(s["S_pct"]): s["ccn_cm3"] for s in post["ccn"]["resolved"][-1]["spectrum"]},
           "nsize": int(run["grid"]["nsize"]), "nclass": int(run["grid"]["nfrac"]),
           "class_labels": run["grid"].get("class_labels")}
    if full:
        out["eabs"] = post["optics"]["resolved"][-1]["bc_abs_enhancement"]
        rr = post["hetchem"]["states"][-1]["rules"]["N2O5"]["resolved"]
        out["k_bt"] = (rr.get("bertram_thornton") or {}).get("k_s1"); out["k_cesm"] = rr["cesm"]["k_s1"]
    out["act"] = {str(w): {**parcel_on(d, w), **fn_arg(d, w)} for w in W}
    return out


def main():
    fo = OUT / f"{sid}.json"
    res = json.loads(fo.read_text()) if fo.exists() else {"id": sid, "configs": {}}
    if "truth" not in res:
        tr = json.loads((SRC / "partmc_s1" / "truth.json").read_text())
        s2 = json.loads((SRC / "stage2.json").read_text())
        res["truth"] = {"ccn": {k: v[-1] for k, v in tr["ccn_cm3"].items()},
                        "nd": {w: s2["activation"][w]["truth"]["nd"] for w in s2["activation"]}}
    for cfg in cfgs:
        if cfg in res["configs"]:
            continue
        t0 = time.time()
        mixing, ov = OC.mixing_and_overrides(cfg)
        d = SCR / cfg.replace("@", "_")
        try:
            PM.project(SRC / "partmc_s1", d, mixing, ov)
            m = metrics(d, full=cfg.endswith("@U1"))
            m["wall_s"] = time.time() - t0
        except Exception as e:          # keep going; record the failure
            m = {"error": repr(e)[:400]}
        res["configs"][cfg] = m
        tmp = fo.with_suffix(".tmp"); tmp.write_text(json.dumps(res)); os.replace(tmp, fo)
        shutil.rmtree(d, ignore_errors=True)
        print(sid, cfg, m.get("nclass"), m.get("nsize"), round(m.get("wall_s", -1), 1), m.get("error", ""), flush=True)
    shutil.rmtree(SCR, ignore_errors=True)


if __name__ == "__main__":
    main()
