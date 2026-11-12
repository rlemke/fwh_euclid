"""Offline tests: tile geometry, the counting query, HiPS layout and the build."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy_healpix import HEALPix
from PIL import Image

from euclidsky import build as bd
from euclidsky import survey as sv
from euclidsky._tap import QueryError, num
from euclidsky.quasars import _sep_arcsec


def _row(t, ra, dec, rel="DR1_R2", first="2024-03-01T00:00:00", last="2024-03-02T00:00:00"):
    return {
        "tile_index": str(t),
        "ra": str(ra),
        "dec": str(dec),
        "release_name": rel,
        "first_observation_date_time": first,
        "last_observation_date_time": last,
        "fov": "(0, 0, 1, 0, 1, 1, 0, 1)",
    }


def test_ring_spacing_is_the_common_gap_not_a_seam_or_a_hole():
    rows = [_row(i, 10 + 0.8 * i, -51.0) for i in range(10)]
    rows.append(_row(99, 10 + 0.8 * 9 + 0.79, -51.0))  # a seam: slightly closer
    rows.append(_row(98, 30.0, -51.0))  # beyond a gap
    tiles = {t.tile_index: t for t in sv.plan_tiles(rows)}
    assert tiles[3].dra == pytest.approx(0.8)
    assert tiles[3].ra_lo == pytest.approx(10 + 2.4 - 0.4)
    assert tiles[3].dec_lo == -51.25


def test_a_reprocessed_tile_keeps_its_later_release():
    rows = [_row(1, 10, 0, "DR1_R1"), _row(1, 10, 0, "DR1_R2")]
    (t,) = sv.plan_tiles(rows)
    assert t.release == "DR1_R2"


def test_core_cut_wraps_through_ra_zero():
    t = sv.Tile(1, 0.1, -60, 359.7, 0.8, -60.25, 0.5, "R", "", "", [])
    where = sv._core_where(t)
    assert "OR" in where and "< 0.500000" in where


def test_counting_query_bins_each_tile_on_its_own_grid():
    t = sv.Tile(7, 80.18, -51, 79.79, 0.786, -51.25, 0.5, "R", "", "", [])
    q = sv.counts_query(t)
    assert "m.tile_index = 7" in q
    assert "/ 0.786000 * 10" in q and "GROUP BY cx, cy, pl, vd, mb" in q
    assert "spurious_flag = 0" in q


def test_nan_and_empty_cells_read_as_missing():
    assert num("") is None and num("NaN") is None and num("1.5") == 1.5


def test_core_area_matches_the_spherical_strip():
    t = {"dra": 0.786, "dec_lo": -51.25, "ddec": 0.5}
    approx = 0.786 * math.cos(math.radians(51)) * 0.5
    assert bd.core_area_deg2(t) == pytest.approx(approx, rel=1e-3)


def test_separation_in_arcsec():
    assert _sep_arcsec(10, 0, 10, 1 / 3600) == pytest.approx(1.0, rel=1e-6)
    assert _sep_arcsec(10, 60, 10 + 2 / 3600, 60) == pytest.approx(1.0, rel=1e-3)


def _tile_file(tmp: Path, *, n: int = 100) -> Path:
    t = sv.Tile(
        5,
        30.0,
        -62.0,
        29.47,
        1.06,
        -62.25,
        0.5,
        "R",
        "2024-01-01",
        "2024-01-02",
        [29.4, -62.27, 30.6, -62.27, 30.6, -61.73, 29.4, -61.73],
    )
    cells = [[cx, cy, 0, 1, 22, n] for cx in range(10) for cy in range(10)]
    d = {"tile": t.__dict__, "grid": 10, "cells": cells, "sample": []}
    (tmp / "tiles").mkdir(parents=True, exist_ok=True)
    f = tmp / "tiles" / "5.json"
    f.write_text(json.dumps(d))
    return f


def test_tile_values_are_densities_per_square_degree(tmp_path):
    d = json.loads(_tile_file(tmp_path).read_text())
    v = bd.tile_values(d)
    area = bd.core_area_deg2(d["tile"]) / 100
    assert v["galaxies"][4, 4] == pytest.approx(100 / area, rel=0.03)
    assert np.isnan(v["stars"]).sum() == 0 and v["stars"].max() == 0


def test_hips_pixel_lands_where_its_sky_position_is(tmp_path):
    """The tile layout Aladin Lite was checked against: a value written for a
    sky position reads back from the PNG at the place the layout says."""
    f = _tile_file(tmp_path)
    vals, _ = bd.pixel_values([f])
    L = bd.LAYERS[0]
    bd.write_hips(tmp_path / "h", L, vals[L.key])
    hp = HEALPix(nside=2**bd.PIX_ORDER, order="nested")
    ip = int(hp.lonlat_to_healpix(30.0 * u.deg, -62.0 * u.deg))
    w = bd.TILE_WIDTH_LOG2
    tile, sub = ip >> (2 * w), ip & ((1 << (2 * w)) - 1)
    px, py = bd._deinterleave(np.array([sub]), w)
    img = np.array(
        Image.open(tmp_path / "h" / f"Norder{bd.TILE_ORDER}" / "Dir0" / f"Npix{tile}.png")
    )
    assert img[(2**w - 1) - px[0], py[0], 3] == 255  # opaque: data here
    ip2 = int(hp.lonlat_to_healpix(30.0 * u.deg, -70.0 * u.deg))
    t2, s2 = ip2 >> (2 * w), ip2 & ((1 << (2 * w)) - 1)
    p2 = tmp_path / "h" / f"Norder{bd.TILE_ORDER}" / "Dir0" / f"Npix{t2}.png"
    if p2.exists():
        x2, y2 = bd._deinterleave(np.array([s2]), w)
        assert np.array(Image.open(p2))[(2**w - 1) - x2[0], y2[0], 3] == 0


def test_build_writes_layers_page_and_summary(tmp_path):
    _tile_file(tmp_path)
    res = bd.build(tmp_path, tmp_path / "site")
    assert res["totals"]["objects"] == 100 * 100
    html = (tmp_path / "site" / "index.html").read_text()
    assert "/*__SUMMARY__*/" not in html and '"objects":10000' in html
    for L in bd.LAYERS:
        assert (tmp_path / "site" / "hips" / L.key / "properties").exists()


def test_query_errors_are_not_retried(monkeypatch):
    from euclidsky import _tap

    calls = []

    def fake(adql, timeout):
        calls.append(1)
        raise QueryError("Unknown column")

    monkeypatch.setattr(_tap, "_once", fake)
    with pytest.raises(QueryError):
        _tap.query("SELECT 1")
    assert len(calls) == 1
