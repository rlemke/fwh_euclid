"""Survey the Euclid DR1 wide catalogue tile by tile, without downloading it.

The archive does the counting. Each MER tile is asked for its objects binned
on a 10 x 10 grid laid over the tile's CORE (tiles overlap their neighbours by
a small margin; counting only the core counts each object once), split by
point-like vs extended, VIS-detected vs NIR-only, and 1-mag brightness bins.
A second, small query takes the tile's brightest extended objects, so the map
can show individual galaxies.

Cells are aligned with the tile, so each is exactly 1/100 of the core in RA
and Dec and its area is known exactly; nothing is resampled at this stage.

Stage outputs, one JSON file per tile, are written atomically, so a run can be
stopped and resumed: a tile whose file exists is not asked again.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import os
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from pathlib import Path

from ._tap import QueryError, num, query

#: Tiles are laid in rings of constant Dec, 0.5 deg apart.
TILE_DEC_STEP = 0.5
GRID = 10
CATALOGUE = "dr1.mer_catalogue_wide_survey"
PHOTO_Z = "dr1.phz_photo_z_wide"


@dataclass
class Tile:
    tile_index: int
    ra: float
    dec: float
    ra_lo: float  # core edge, [0, 360)
    dra: float  # core width in RA, degrees of RA
    dec_lo: float
    ddec: float
    release: str
    first_obs: str
    last_obs: str
    fov: list[float]  # footprint polygon as published (ra, dec pairs)


def plan_tiles(rows: list[dict[str, str]]) -> list[Tile]:
    """Wide-survey tiles from ``dr1.mosaic_product`` VIS rows, with each tile's
    core. The RA width of a ring's tiles is the spacing between neighbours in
    that ring; a tile with no neighbour in its ring gets the nominal width."""
    by_tile: dict[int, dict[str, str]] = {}
    for r in rows:
        t = int(r["tile_index"])
        # DR1_R2 reprocessed some R1 tiles; keep the later one.
        if t not in by_tile or r["release_name"] > by_tile[t]["release_name"]:
            by_tile[t] = r
    rings: dict[str, list[float]] = {}
    for r in by_tile.values():
        rings.setdefault(r["dec"], []).append(float(r["ra"]))
    spacing: dict[str, float] = {}
    for d, ras in rings.items():
        ras = sorted(set(ras))
        diffs = [round(b - a, 4) for a, b in itertools.pairwise(ras) if b - a > 1e-6]
        nominal = TILE_DEC_STEP / max(math.cos(math.radians(float(d))), 0.05)
        # The ring's own spacing is its most common neighbour gap: seams
        # between survey patches give odd gaps, and gaps over unobserved sky
        # give multiples.
        near = [x for x in diffs if x < 1.5 * nominal]
        spacing[d] = Counter(near).most_common(1)[0][0] if near else nominal
    out = []
    for t, r in sorted(by_tile.items()):
        ra, dec = float(r["ra"]), float(r["dec"])
        dra = spacing[r["dec"]]
        fov = [float(x) for x in r["fov"].strip("()").split(",")]
        out.append(
            Tile(
                tile_index=t,
                ra=ra,
                dec=dec,
                ra_lo=(ra - dra / 2) % 360,
                dra=dra,
                dec_lo=dec - TILE_DEC_STEP / 2,
                ddec=TILE_DEC_STEP,
                release=r["release_name"],
                first_obs=r["first_observation_date_time"],
                last_obs=r["last_observation_date_time"],
                fov=fov,
            )
        )
    return out


def fetch_tile_list() -> list[dict[str, str]]:
    return query(
        "SELECT tile_index, ra, dec, release_name, first_observation_date_time, "
        "last_observation_date_time, fov FROM dr1.mosaic_product "
        "WHERE filter_name = 'VIS' AND processing_mode = 'WIDE'"
    )


def _core_where(t: Tile, a: str = "m") -> str:
    ra_hi = t.ra_lo + t.dra
    dec = f"{a}.declination >= {t.dec_lo:.6f} AND {a}.declination < {t.dec_lo + t.ddec:.6f}"
    if ra_hi <= 360:
        ra = f"{a}.right_ascension >= {t.ra_lo:.6f} AND {a}.right_ascension < {ra_hi:.6f}"
    else:
        ra = f"({a}.right_ascension >= {t.ra_lo:.6f} OR {a}.right_ascension < {ra_hi - 360:.6f})"
    return f"{dec} AND {ra}"


def counts_query(t: Tile) -> str:
    cx = f"FLOOR(MOD(m.right_ascension - {t.ra_lo:.6f} + 720, 360) / {t.dra:.6f} * {GRID})"
    cy = f"FLOOR((m.declination - {t.dec_lo:.6f}) / {t.ddec:.6f} * {GRID})"
    return (
        "SELECT cx, cy, pl, vd, mb, COUNT(*) AS n FROM ("
        f"SELECT {cx} AS cx, {cy} AS cy, FLOOR(m.point_like_prob * 2) AS pl, "
        "m.vis_det AS vd, FLOOR(23.9 - 2.5 * LOG10(m.flux_detection_total)) AS mb "
        f"FROM {CATALOGUE} m WHERE m.tile_index = {t.tile_index} "
        f"AND m.spurious_flag = 0 AND m.flux_detection_total > 0 AND {_core_where(t)}"
        ") AS s GROUP BY cx, cy, pl, vd, mb"
    )


def bright_query(t: Tile, n: int = 150) -> str:
    """The tile's brightest objects whatever their shape: the point-like
    probability is blank for saturated and NIR-only objects, which is where
    the biggest galaxies are, so the class comes from the photometric
    classifier instead (see ``sample_tile``)."""
    return (
        f"SELECT TOP {n} m.object_id, m.right_ascension, m.declination, "
        "m.flux_detection_total, m.semimajor_axis, m.ellipticity, m.position_angle, m.gaia_id, "
        "m.mumax_minus_mag "
        f"FROM {CATALOGUE} m WHERE m.tile_index = {t.tile_index} "
        "AND m.spurious_flag = 0 AND m.flux_detection_total > 0 "
        f"AND m.flux_detection_total < 100000000 AND {_core_where(t)} "
        "ORDER BY m.flux_detection_total DESC"
    )


def sample_tile(t: Tile, keep: int = 30) -> dict:
    """The brightest objects Euclid's photometric classifier calls galaxies
    (class flag exactly 2), with no Gaia counterpart and a diffuse profile."""
    rows = query(bright_query(t))
    out = []
    if rows:
        zs = {r["object_id"]: r for r in query(z_query([r["object_id"] for r in rows]))}
        for r in rows:
            z = zs.get(r["object_id"], {})
            if _int(z.get("phz_classification", "")) != 2 or r.get("gaia_id"):
                continue
            # Saturated stars can pass the classifier; a star's peak surface
            # brightness sits within about a magnitude of its total magnitude,
            # a galaxy's well below it.
            mm = num(r.get("mumax_minus_mag"))
            if mm is None or mm < 1.0:
                continue
            out.append(
                {
                    "ra": float(r["right_ascension"]),
                    "dec": float(r["declination"]),
                    "flux": num(r["flux_detection_total"]),
                    "a_px": num(r["semimajor_axis"]),
                    "e": num(r["ellipticity"]),
                    "pa": num(r["position_angle"]),
                    "z": num(z.get("phz_median")),
                    "cls": 2,
                }
            )
            if len(out) >= keep:
                break
    return {"tile_index": t.tile_index, "sample": out}


def z_query(ids: list[str]) -> str:
    return (
        "SELECT object_id, phz_median, phz_classification "
        f"FROM {PHOTO_Z} WHERE object_id IN ({', '.join(ids)})"
    )


def _int(s: str) -> int | None:
    v = num(s)
    return None if v is None else int(v)


def survey_tile(t: Tile) -> dict:
    """One tile's counts and brightest galaxies, as a JSON-able dict."""
    t0 = time.time()
    rows = query(counts_query(t))
    cells = []
    for r in rows:
        cx, cy = _int(r["cx"]), _int(r["cy"])
        if cx is None or cy is None:
            continue
        cells.append(
            [
                min(max(cx, 0), GRID - 1),
                min(max(cy, 0), GRID - 1),
                _int(r["pl"]),  # 0 extended, 1-2 point-like, None no VIS shape
                _int(r["vd"]),
                _int(r["mb"]),
                int(r["n"]),
            ]
        )
    # The galaxy sample is its own stage (``sample_tile``).
    sample: list[dict] = []
    return {
        "tile": asdict(t),
        "grid": GRID,
        "cells": cells,
        "sample": sample,
        "seconds": round(time.time() - t0, 1),
    }


def _write_atomic(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(obj, separators=(",", ":")))
    os.replace(tmp, path)


def load_tiles(out: Path) -> list[Tile]:
    """The planned tiles, fetching the archive's list once into ``out``."""
    out.mkdir(parents=True, exist_ok=True)
    tl = out / "tiles.csv"
    if not tl.exists():
        rows = fetch_tile_list()
        tmp = tl.with_suffix(".tmp")
        with tmp.open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)
        os.replace(tmp, tl)
    return plan_tiles(list(csv.DictReader(tl.open())))


def run(
    out: Path,
    workers: int,
    limit: int | None = None,
    stage: str = "counts",
    only: set[int] | None = None,
) -> int:
    """Run a stage over every planned tile (or ``only`` these), skipping
    tiles whose file exists. Returns 1 if any tile failed."""
    tiles = load_tiles(out)
    if only is not None:
        tiles = [t for t in tiles if t.tile_index in only]
    sub = "tiles" if stage == "counts" else "samples"
    fn = survey_tile if stage == "counts" else sample_tile
    (out / sub).mkdir(exist_ok=True)
    if stage == "samples":
        # Only tiles the counts stage found objects in.
        have = set()
        for p in (out / "tiles").glob("*.json"):
            if json.loads(p.read_text())["cells"]:
                have.add(int(p.stem))
        tiles = [t for t in tiles if t.tile_index in have]
    todo = [t for t in tiles if not (out / sub / f"{t.tile_index}.json").exists()]
    if limit:
        todo = todo[:limit]
    print(f"{len(tiles)} tiles, {len(todo)} to fetch, {workers} workers", flush=True)
    failed = 0
    done = 0
    t0 = time.time()
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(fn, t): t for t in todo}
        for f in as_completed(futs):
            t = futs[f]
            try:
                res = f.result()
            except QueryError as e:
                failed += 1
                print(f"tile {t.tile_index}: rejected: {e}", file=sys.stderr, flush=True)
                continue
            except Exception as e:  # noqa: BLE001 - report and carry on; rerun resumes
                failed += 1
                print(f"tile {t.tile_index}: failed: {e}", file=sys.stderr, flush=True)
                continue
            _write_atomic(out / sub / f"{t.tile_index}.json", res)
            done += 1
            if done % 100 == 0:
                rate = done / (time.time() - t0)
                left = (len(todo) - done) / rate / 3600 if rate else 0
                print(f"{done}/{len(todo)} tiles, {left:.1f} h left", flush=True)
    print(f"done: {done} fetched, {failed} failed", flush=True)
    return 1 if failed else 0


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="euclidsky.survey", description=__doc__.split("\n")[0])
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--limit", type=int)
    p.add_argument("--stage", choices=("counts", "samples"), default="counts")
    a = p.parse_args(argv)
    return run(a.out, a.workers, a.limit, a.stage)


if __name__ == "__main__":
    sys.exit(main())
