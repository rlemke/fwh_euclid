"""How Euclid DR1 labels quasars that other telescopes have already confirmed.

Reference: the Million Quasars catalogue (Milliquas, Flesch), entries whose
type says a SPECTRUM confirmed them (Q = type-1 quasar; A, B, K, N = other
spectroscopically confirmed active nuclei). Photometric candidates in
Milliquas (R, X, 2) are not used as truth.

For every confirmed quasar inside a DR1 tile core, the archive is asked for
the Euclid object within 1 arcsec, and then for each of Euclid's labels on
that object:

* PHZ class flag (bit 4 = accepted as a QSO by the photometric classifier),
  with the QSO probability and the threshold it was compared against;
* whether a QSO template fit exists for it (``phz_qso_physical_parameters``);
* the spectroscopic class from Euclid's own slitless spectrum, and the
  quasar-template spectroscopic redshift.

This measures COMPLETENESS: of the quasars known to be real, how many does
Euclid call quasars, and how well do its redshifts agree. Purity, the other
half, needs every object Euclid calls a quasar to have been checked, which
no reference catalogue does; ``purity`` estimates it where the reference is
deepest, as a lower bound.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from ._tap import num, query

CATALOGUES = ("dr1.mer_catalogue_wide_survey", "dr1.mer_catalogue_wide_mode")
MATCH_ARCSEC = 1.0


def _ids(rows: list[dict], key: str = "object_id") -> str:
    return ", ".join(r[key] for r in rows)


def match_tile(tile: int, refs: list[dict]) -> list[dict]:
    """The Euclid object nearest each reference position (within 1"). A
    reference with no Euclid object is returned unmatched; labels come later
    (``labels``), batched across tiles."""
    r_deg = 1.5 / 3600
    boxes = []
    for r in refs:
        ra, dec = float(r["ra"]), float(r["dec"])
        dra = r_deg / max(math.cos(math.radians(dec)), 0.01)
        boxes.append(
            f"(m.right_ascension BETWEEN {ra - dra:.7f} AND {ra + dra:.7f} "
            f"AND m.declination BETWEEN {dec - r_deg:.7f} AND {dec + r_deg:.7f})"
        )
    cand: list[dict] = []
    cat_used = None
    for cat in CATALOGUES:
        cand = query(
            "SELECT m.object_id, m.right_ascension, m.declination, m.flux_detection_total, "
            "m.flux_vis_psf, m.point_like_prob, m.vis_det, m.spurious_flag "
            f"FROM {cat} m WHERE m.tile_index = {tile} AND ({' OR '.join(boxes)})"
        )
        if cand:
            cat_used = cat
            break
    out = []
    for r in refs:
        ra, dec = float(r["ra"]), float(r["dec"])
        best, best_d = None, MATCH_ARCSEC
        for c in cand:
            d = _sep_arcsec(ra, dec, float(c["right_ascension"]), float(c["declination"]))
            if d <= best_d:
                best, best_d = c, d
        row = {"ref": r, "tile": tile, "catalogue": cat_used, "match": None}
        if best:
            row["match"] = {"sep": round(best_d, 3), **best}
        out.append(row)
    return out


def labels(object_ids: list[str]) -> dict[str, dict]:
    """Every Euclid label for these objects, a few queries for the lot: each
    archive query costs ~2 s whatever its size, so labels are fetched in
    batches across tiles rather than per tile."""
    ids = ", ".join(object_ids)
    lab: dict[str, dict] = defaultdict(dict)
    for r in query(
        "SELECT object_id, phz_classification, phz_qso_prob, qso_threshold, phz_star_prob, "
        f"phz_gal_prob FROM dr1.phz_classification_wide WHERE object_id IN ({ids})"
    ):
        lab[r["object_id"]].update(r)
    for r in query(
        f"SELECT object_id, phz_median FROM dr1.phz_photo_z_wide WHERE object_id IN ({ids})"
    ):
        lab[r["object_id"]]["phz_median"] = r["phz_median"]
    for r in query(
        f"SELECT object_id, z FROM dr1.phz_qso_physical_parameters_wide WHERE object_id IN ({ids})"
    ):
        lab[r["object_id"]]["qso_fit_z"] = r["z"]
    for r in query(
        "SELECT object_id, spe_class, spe_qso_prob FROM "
        f"dr1.spectro_z_spe_classification_wide_survey WHERE object_id IN ({ids})"
    ):
        lab[r["object_id"]].update(r)
    for r in query(
        "SELECT object_id, spe_z, spe_z_prob FROM dr1.spectro_z_spe_qso_candidates_wide_survey "
        f"WHERE spe_rank = 0 AND object_id IN ({ids})"
    ):
        lab[r["object_id"]]["spe_qso_z"] = r["spe_z"]
    lab_galaxy_z(lab, ids)
    for i in object_ids:
        lab[i]["labelled"] = "1"
    return dict(lab)


def lab_galaxy_z(lab: dict, ids: str) -> None:
    """The galaxy-template spectroscopic redshift: what Euclid reports for an
    object its spectrum calls a galaxy."""
    for r in query(
        "SELECT object_id, spe_z FROM dr1.spectro_z_spe_galaxy_candidates_wide_survey "
        f"WHERE spe_rank = 0 AND object_id IN ({ids})"
    ):
        lab[r["object_id"]]["spe_gal_z"] = r["spe_z"]


def spectral_qsos(tile: int) -> list[dict]:
    """Every object in a tile that Euclid's spectrum classifies as a quasar,
    with position and brightness, for the reverse check against the
    reference."""
    return query(
        "SELECT m.object_id, m.right_ascension, m.declination, m.flux_detection_total, "
        "m.point_like_prob, s.spe_qso_prob FROM dr1.mer_catalogue_wide_survey m "
        "JOIN dr1.spectro_z_spe_classification_wide_survey s ON s.object_id = m.object_id "
        f"WHERE m.tile_index = {tile} AND s.spe_class = 'qso'"
    )


def _sep_arcsec(ra1: float, d1: float, ra2: float, d2: float) -> float:
    r1, r2 = math.radians(d1), math.radians(d2)
    dra = math.radians(ra2 - ra1)
    a = math.sin((r2 - r1) / 2) ** 2 + math.cos(r1) * math.cos(r2) * math.sin(dra / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(a)))) * 3600


def run_reverse(refs_csv: Path, out: Path, workers: int, n_tiles: int = 200) -> int:
    """Euclid's spectroscopic quasars in the tiles where the reference is
    densest, so a missing reference match says most about Euclid rather than
    about the reference."""
    per: dict[int, int] = defaultdict(int)
    for r in csv.DictReader(refs_csv.open()):
        if r["type"][:1] == "Q":
            per[int(r["tile"])] += 1
    done = {p.stem for p in (out / "xmatch").glob("*.json")}
    survey_tiles = []
    for t, n in sorted(per.items(), key=lambda x: -x[1]):
        f = out / "xmatch" / f"{t}.json"
        if str(t) in done and json.loads(f.read_text())[0]["catalogue"] == CATALOGUES[0]:
            survey_tiles.append(t)
        if len(survey_tiles) >= n_tiles:
            break
    (out / "reverse").mkdir(parents=True, exist_ok=True)
    todo = [t for t in survey_tiles if not (out / "reverse" / f"{t}.json").exists()]
    print(f"reverse check over {len(survey_tiles)} tiles, {len(todo)} to do", flush=True)
    failed = 0
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(spectral_qsos, t): t for t in todo}
        for f in as_completed(futs):
            t = futs[f]
            try:
                rows = f.result()
            except Exception as e:  # noqa: BLE001 - rerun resumes
                failed += 1
                print(f"tile {t}: {str(e)[:300]}", file=sys.stderr, flush=True)
                continue
            (out / "reverse" / f"{t}.json").write_text(json.dumps({"tile": t, "rows": rows}))
    print(f"reverse done, {failed} failed", flush=True)
    return 1 if failed else 0


def run(refs_csv: Path, out: Path, workers: int) -> int:
    rows = [
        r for r in csv.DictReader(refs_csv.open()) if r["type"][:1] in ("Q", "A", "B", "K", "N")
    ]
    by_tile: dict[int, list[dict]] = defaultdict(list)
    for r in rows:
        by_tile[int(r["tile"])].append(r)
    (out / "xmatch").mkdir(parents=True, exist_ok=True)
    todo = [t for t in by_tile if not (out / "xmatch" / f"{t}.json").exists()]
    print(f"{len(rows)} confirmed AGN in {len(by_tile)} tiles; {len(todo)} to do", flush=True)
    failed = 0
    t0 = time.time()
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(match_tile, t, by_tile[t]): t for t in todo}
        for i, f in enumerate(as_completed(futs), 1):
            t = futs[f]
            try:
                res = f.result()
            except Exception as e:  # noqa: BLE001 - rerun resumes
                failed += 1
                print(f"tile {t}: {str(e)[:300]}", file=sys.stderr, flush=True)
                continue
            p = out / "xmatch" / f"{t}.json"
            tmp = p.with_suffix(".tmp")
            tmp.write_text(json.dumps(res))
            os.replace(tmp, p)
            if i % 100 == 0:
                print(f"{i}/{len(todo)} tiles, {(time.time() - t0) / 60:.0f} min", flush=True)
    print(f"matched, {failed} failed; fetching labels", flush=True)
    failed += run_labels(out, workers)
    print(f"done, {failed} failed", flush=True)
    return 1 if failed else 0


def run_labels(out: Path, workers: int, batch: int = 400) -> int:
    """Labels for every matched object not yet labelled, in batches."""
    have: set[str] = set()
    ldir = out / "labels"
    ldir.mkdir(exist_ok=True)
    for p in ldir.glob("*.json"):
        have.update(json.loads(p.read_text()))
    ids = []
    for p in sorted((out / "xmatch").glob("*.json")):
        for o in json.loads(p.read_text()):
            m = o["match"]
            if m and o["catalogue"] == CATALOGUES[0] and m["object_id"] not in have:
                ids.append(m["object_id"])
    ids = sorted(set(ids))
    batches = [ids[i : i + batch] for i in range(0, len(ids), batch)]
    print(f"{len(ids)} objects to label in {len(batches)} batches", flush=True)
    failed = 0
    start = len(list(ldir.glob("*.json")))
    with ThreadPoolExecutor(workers) as ex:
        futs = {ex.submit(labels, b): i for i, b in enumerate(batches)}
        for f in as_completed(futs):
            i = futs[f]
            try:
                lab = f.result()
            except Exception as e:  # noqa: BLE001 - rerun resumes
                failed += 1
                print(f"label batch {i}: {str(e)[:300]}", file=sys.stderr, flush=True)
                continue
            (ldir / f"{start + i:05d}.json").write_text(json.dumps(lab))
    return failed


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="euclidsky.quasars")
    p.add_argument("--refs", type=Path, required=True, help="reference quasars with a tile column")
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--stage", choices=("match", "reverse"), default="match")
    a = p.parse_args(argv)
    if a.stage == "reverse":
        return run_reverse(a.refs, a.out, a.workers)
    return run(a.refs, a.out, a.workers)


if __name__ == "__main__":
    sys.exit(main())


__all__ = ["match_tile", "num"]
