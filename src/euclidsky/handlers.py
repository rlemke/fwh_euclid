"""Facetwork handlers for the ``euclid`` facets.

Each handler is a thin step over the code the command-line tools run, so a
workflow run and a CLI run write the same files.

⚠️ ``work_dir`` and ``site_dir`` are LOCAL paths: the per-tile files are
small and many, and the build reads all of them. A fleet run needs every
step on one host (or a shared volume), which is why the domain has no
``service`` and runs from a checkout.
"""

from __future__ import annotations

import os
from collections.abc import Callable
from pathlib import Path
from typing import Any

from . import build as bd
from . import survey as sv


def _say(params: dict[str, Any]) -> Callable[..., None]:
    cb = params.get("_step_log")
    return cb if callable(cb) else (lambda *a, **k: None)


def handle_plan(params: dict[str, Any]) -> dict[str, Any]:
    out = Path(params["work_dir"])
    size = int(params.get("chunk_size") or 250)
    tiles = sv.load_tiles(out)
    ids = [t.tile_index for t in tiles]
    chunks = [ids[i : i + size] for i in range(0, len(ids), size)]
    _say(params)(f"{len(ids)} tiles in {len(chunks)} chunks of up to {size}")
    return {"chunks": chunks, "tiles": len(ids)}


def _stage(params: dict[str, Any], stage: str) -> tuple[int, int]:
    out = Path(params["work_dir"])
    only = {int(t) for t in params["tiles"]}
    sv.run(out, int(params.get("workers") or 4), stage=stage, only=only)
    sub = "tiles" if stage == "counts" else "samples"
    have = {int(p.stem) for p in (out / sub).glob("*.json")}
    done = len(only & have)
    if stage == "samples":
        # Only tiles with objects are sampled.
        with_objects = 0
        for t in only:
            p = out / "tiles" / f"{t}.json"
            if p.exists() and '"cells":[]' not in p.read_text()[:4000]:
                with_objects += 1
        return done, max(0, with_objects - done)
    return done, len(only) - done


def handle_count_chunk(params: dict[str, Any]) -> dict[str, Any]:
    counted, failed = _stage(params, "counts")
    _say(params)(f"counted {counted} tiles, {failed} failed (a retry resumes them)")
    return {"counted": counted, "failed": failed}


def handle_sample_chunk(params: dict[str, Any]) -> dict[str, Any]:
    sampled, failed = _stage(params, "samples")
    _say(params)(f"sampled {sampled} tiles, {failed} failed")
    return {"sampled": sampled, "failed": failed}


def handle_build_sky(params: dict[str, Any]) -> dict[str, Any]:
    res = bd.build(Path(params["work_dir"]), Path(params["site_dir"]))
    t = res["totals"]
    _say(params)(f"{t['objects']:,} objects over {t['area_deg2']:,.0f} deg2")
    return {
        "index": str(Path(params["site_dir"]) / "index.html"),
        "objects": t["objects"],
        "galaxies": t["galaxies"],
        "area_deg2": t["area_deg2"],
    }


_DISPATCH: dict[str, Callable[[dict[str, Any]], dict[str, Any]]] = {
    "euclid.Plan": handle_plan,
    "euclid.CountChunk": handle_count_chunk,
    "euclid.SampleChunk": handle_sample_chunk,
    "euclid.BuildSky": handle_build_sky,
}


def handle(payload: dict[str, Any]) -> dict[str, Any]:
    """RegistryRunner entrypoint: route on the facet name the runner supplies."""
    fn = _DISPATCH.get(payload.get("_facet_name", ""))
    if fn is None:
        raise ValueError(f"no euclid handler for facet {payload.get('_facet_name')!r}")
    return fn(payload)


def register_handlers(runner: Any) -> None:
    """RegistryRunner registration -- one entrypoint per facet."""
    for facet_name in _DISPATCH:
        runner.register_handler(
            facet_name=facet_name,
            module_uri=f"file://{os.path.abspath(__file__)}",
            entrypoint="handle",
        )


def register_poller(poller: Any) -> None:
    for facet_name, fn in _DISPATCH.items():
        poller.register(facet_name, fn)
