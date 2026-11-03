"""Post-processing identical to SCRAM: optics/radiation, cloud activation (CCN),
heterogeneous uptake. Each module works on one output state (common.Populations) and
reports the class-resolved result next to the internally mixed reference."""
from __future__ import annotations

import json
from pathlib import Path

from . import ccn, hetchem, optics

MODULES = {"optics": optics, "ccn": ccn, "hetchem": hetchem}


def run_all(run: dict, cfg: dict | None = None, run_dir: str | Path | None = None) -> dict:
    """Run the enabled modules (cfg = case.postprocess); write post.json into run_dir."""
    cfg = cfg if cfg is not None else run["case"].get("postprocess", {})
    out = {}
    for name, mod in MODULES.items():
        c = dict(cfg.get(name) or {})
        if not c.pop("enabled", False):
            continue
        every = int(c.pop("every", 1))
        out[name] = mod.series(run, c, every)
    if run_dir is not None and out:
        Path(run_dir, "post.json").write_text(json.dumps(out, ensure_ascii=False))
    return out
