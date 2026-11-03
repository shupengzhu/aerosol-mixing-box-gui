"""Command line: python -m scrambox <command> ...

  run CASE.json -o RUN_DIR      compile and run a case
  compile CASE.json -o RUN_DIR  compile only (inspect run.nml / inputs.txt)
  classes [--preset NC6]        list the composition classes of the default species table
  serve [--port 8765]           start the web interface
"""
from __future__ import annotations

import argparse
import json
import sys

from . import case as casemod


def main(argv=None):
    ap = argparse.ArgumentParser(prog="scrambox")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "compile"):
        p = sub.add_parser(name)
        p.add_argument("case")
        p.add_argument("-o", "--out", required=True)
    p = sub.add_parser("classes")
    p.add_argument("--preset", default="NC6")
    p.add_argument("--isorropia", default="v23")
    p = sub.add_parser("serve")
    p.add_argument("--port", type=int, default=8765)
    p.add_argument("--host", default="127.0.0.1")
    a = ap.parse_args(argv)

    if a.cmd == "compile":
        from .compile import compile_case
        info = compile_case(casemod.load(a.case), a.out)
        print(json.dumps({"run_dir": info["run_dir"], "nsize": info["grid"]["nsize"],
                          "nfrac": info["grid"]["nfrac"], "notes": info["grid"]["notes"]}, indent=1))
    elif a.cmd == "run":
        from .run import run_case
        out = run_case(casemod.load(a.case), a.out)
        print(json.dumps(out["summary"], indent=1, ensure_ascii=False))
    elif a.cmd == "classes":
        from .compile import build_grid
        c = casemod.resolve({"mixing": {"preset": a.preset}, "kernel": {"isorropia": a.isorropia}})
        g = build_grid(c)
        from .grid import class_label
        print("groups:", g["table"].groups)
        for i, cls in enumerate(g["classes"], 1):
            print(f"C{i:<3d}", class_label(cls, g["table"].groups), cls)
    elif a.cmd == "serve":
        from .server import serve
        serve(a.host, a.port)


if __name__ == "__main__":
    sys.exit(main())
