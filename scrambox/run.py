"""Run a compiled case with the matching kernel executable and read the results."""
from __future__ import annotations

import json
import os
import re
import subprocess
import time
from pathlib import Path

import numpy as np

from .compile import compile_case

ROOT = Path(__file__).resolve().parent.parent
KERNEL_BUILD = Path(os.environ.get("SCRAMBOX_KERNEL_BUILD", ROOT / "kernel" / "build"))


def executable(isorropia: str) -> Path:
    """v23 kernel for ISORROPIA II 2.3; v1 kernel for ISORROPIA v1 and for no-ISORROPIA."""
    exe = KERNEL_BUILD / ("v23" if isorropia == "v23" else "v1") / "scram_box"
    if not exe.exists():
        raise FileNotFoundError(f"kernel executable not found: {exe} (run kernel/build_kernel.sh)")
    return exe


def kernel_version(isorropia: str) -> str:
    f = KERNEL_BUILD / ("v23" if isorropia == "v23" else "v1") / "KERNEL_VERSION"
    return f.read_text() if f.exists() else "unknown"


_FEXP = re.compile(r"(?<=\d)([+-]\d{3})$")


def _num(tok: str) -> float:
    """Fortran ES output drops the 'E' for 3-digit exponents (1.0-112); accept both."""
    try:
        return float(tok)
    except ValueError:
        return float(_FEXP.sub(r"E\1", tok))


def read_output(path: str | Path, nsize: int, nfrac: int, naer: int, ngas: int) -> dict:
    nbin = nsize * nfrac
    times, met, num, mass, gas = [], [], [], [], []
    cpu = None
    with open(path) as f:
        tokens = None
        for line in f:
            parts = line.split()
            if not parts or parts[0].startswith("#"):
                continue
            key = parts[0]
            if key == "time":
                times.append(_num(parts[1]))
                met.append([_num(x) for x in parts[2:5]])
            elif key in ("number_m3", "mass_ugm3", "gas_ugm3"):
                n = int(parts[1])
                vals = []
                while len(vals) < n:
                    vals += [_num(x) for x in next(f).split()]
                a = np.array(vals)
                if key == "number_m3":
                    num.append(a)
                elif key == "mass_ugm3":
                    mass.append(a.reshape((nbin, naer), order="F"))
                else:
                    gas.append(a)
            elif key == "cpu_s":
                cpu = _num(parts[1])
    num = np.array(num).reshape(len(times), nsize, nfrac)
    mass = np.array(mass).reshape(len(times), nsize, nfrac, naer)
    return {"time_s": np.array(times), "met": np.array(met), "number_m3": num,
            "mass_ugm3": mass, "gas_ugm3": np.array(gas), "cpu_s": cpu}


def conservation(res: dict, grid: dict) -> dict:
    """Relative change of gas + particle mass for each condensable pair, first to last
    output (meaningful when emissions and gas sources are zero)."""
    sp = grid["species"]
    out = {}
    for aer, gasname in sp["gas_of"].items():
        if aer not in sp["aerosol"] or gasname not in sp["gases"]:
            continue
        ia, ig = sp["aerosol"].index(aer), sp["gases"].index(gasname)
        tot = res["mass_ugm3"][:, :, :, ia].sum(axis=(1, 2)) + res["gas_ugm3"][:, ig]
        if tot[0] > 0:
            out[f"{aer}+{gasname}"] = float((tot[-1] - tot[0]) / tot[0])
    return out


def run_case(user_case: dict, run_dir: str | Path, timeout_s: float = 7200) -> dict:
    run = Path(run_dir).resolve()
    info = compile_case(user_case, run)
    c, g = info["case"], info["grid"]
    env = dict(os.environ)
    env["SCRAM_SPECIES_FILE"] = str(run / "scram_species.dat")
    env.setdefault("OMP_NUM_THREADS", "1")
    t0 = time.time()
    with open(run / "run.log", "w") as log:
        proc = subprocess.run([str(executable(c["kernel"]["isorropia"])), "run.nml"], cwd=run,
                              env=env, stdout=log, stderr=subprocess.STDOUT, timeout=timeout_s)
    wall = time.time() - t0
    summary = {"name": c["case"]["name"], "returncode": proc.returncode, "wall_s": wall,
               "kernel": kernel_version(c["kernel"]["isorropia"]), "notes": g["notes"]}
    if proc.returncode != 0:
        summary["error"] = (run / "run.log").read_text()[-2000:]
        (run / "summary.json").write_text(json.dumps(summary, indent=1))
        raise RuntimeError(f"scram_box failed (rc={proc.returncode}): {summary['error'][-400:]}")
    sp = g["species"]
    res = read_output(run / "output.txt", g["nsize"], g["nfrac"], len(sp["aerosol"]), len(sp["gases"]))
    np.savez_compressed(run / "results.npz", **{k: v for k, v in res.items() if v is not None})
    summary.update({"cpu_s": res["cpu_s"], "n_output": int(len(res["time_s"])),
                    "conservation_rel": conservation(res, g)})
    out = {"case": c, "grid": g, "results": res, "summary": summary}
    if any((c.get("postprocess") or {}).get(k, {}).get("enabled") for k in ("optics", "ccn", "hetchem")):
        from .post import run_all
        t1 = time.time()
        out["post"] = run_all(out, c["postprocess"], run)
        summary["post_wall_s"] = time.time() - t1
    (run / "summary.json").write_text(json.dumps(summary, indent=1, ensure_ascii=False))
    return out


def load_run(run_dir: str | Path) -> dict:
    run = Path(run_dir)
    g = json.loads((run / "grid.json").read_text())
    data = np.load(run / "results.npz")
    res = {k: data[k] for k in data.files}
    out = {"case": json.loads((run / "case.resolved.json").read_text()), "grid": g,
           "results": res, "summary": json.loads((run / "summary.json").read_text())}
    if (run / "post.json").exists():
        out["post"] = json.loads((run / "post.json").read_text())
    return out
