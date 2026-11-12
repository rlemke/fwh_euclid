"""Turn the per-tile survey files into sky layers and a summary.

Layers are written as HiPS (the IVOA tiled all-sky format sky viewers read),
so the density shows on the sky itself, over real imagery, at the resolution
the counts were taken: a 10 x 10 grid per tile core, about 3 arcmin.

Each HiPS pixel takes the value of the survey cell its centre falls in. That
is a lookup, not a re-binning: the counts were never moved between cells, so
no cell's density is spread into its neighbours or aliased onto the pixel
grid. Pixels outside every tile core are transparent.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy_healpix import HEALPix
from PIL import Image

#: Finest HiPS pixel order (1.7 arcmin), below the ~3 arcmin survey cells.
PIX_ORDER = 11
TILE_WIDTH_LOG2 = 9  # 512 x 512 pixel tiles: Aladin Lite draws nothing for 128
TILE_ORDER = PIX_ORDER - TILE_WIDTH_LOG2  # 2

#: The uniform-depth galaxy cut: extended objects brighter than this
#: (detection-band total magnitude). Euclid wide reaches ~24.5 in VIS at 10
#: sigma, so 24 is complete over the whole footprint where depth is nominal.
GAL_MAG_LIMIT = 24


@dataclass
class Layer:
    key: str
    title: str
    unit: str
    vmin: float
    vmax: float
    log: bool
    cmap: str


SAMPLE_KEEP = 30


def core_area_deg2(t: dict) -> float:
    d0, d1 = math.radians(t["dec_lo"]), math.radians(t["dec_lo"] + t["ddec"])
    return math.radians(t["dra"]) * (math.sin(d1) - math.sin(d0)) * (180 / math.pi) ** 2


def _cell_area_deg2(t: dict, cy: int, grid: int) -> float:
    d0 = math.radians(t["dec_lo"] + cy * t["ddec"] / grid)
    d1 = math.radians(t["dec_lo"] + (cy + 1) * t["ddec"] / grid)
    dra = math.radians(t["dra"] / grid)
    return dra * (math.sin(d1) - math.sin(d0)) * (180 / math.pi) ** 2


def tile_values(d: dict) -> dict[str, np.ndarray]:
    """Per-cell values (grid x grid arrays, NaN = no data) for every layer."""
    t, g = d["tile"], d["grid"]
    allobj = np.zeros((g, g))
    gal = np.zeros((g, g))
    star = np.zeros((g, g))
    for cx, cy, pl, vd, mb, n in d["cells"]:
        allobj[cy, cx] += n
        if pl == 0 and mb is not None and mb < GAL_MAG_LIMIT:
            gal[cy, cx] += n
        if pl is not None and pl >= 1:
            star[cy, cx] += n
    area = np.array([[_cell_area_deg2(t, cy, g)] * g for cy in range(g)])
    empty = allobj == 0
    out = {}
    for k, a in (("all", allobj), ("galaxies", gal), ("stars", star)):
        v = a / area
        v[empty] = np.nan
        out[k] = v
    obs = _obs_days(t.get("first_obs"))
    out["observed"] = np.where(empty, np.nan, obs if obs is not None else np.nan)
    return out


SURVEY_START = date(2023, 10, 1)


def _obs_days(s: str | None) -> float | None:
    if not s:
        return None
    try:
        return (datetime.fromisoformat(s[:19]).date() - SURVEY_START).days
    except ValueError:
        return None


LAYERS = [
    Layer(
        "galaxies",
        f"Galaxies brighter than mag {GAL_MAG_LIMIT}",
        "per deg²",
        2e4,
        2.2e5,
        True,
        "magma",
    ),
    Layer("all", "All detected objects", "per deg²", 5e4, 6e5, True, "viridis"),
    Layer(
        "stars",
        "Point-like objects (mostly Milky Way stars)",
        "per deg²",
        1e3,
        2e5,
        True,
        "cividis",
    ),
    Layer(
        "observed",
        "When Euclid first observed the tile",
        "days since 2023-10-01",
        0,
        600,
        False,
        "turbo",
    ),
]


def _cmap(name: str):
    from matplotlib import colormaps

    return colormaps[name]


def _deinterleave(s: np.ndarray, bits: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.zeros_like(s)
    y = np.zeros_like(s)
    for b in range(bits):
        x |= ((s >> (2 * b)) & 1) << b
        y |= ((s >> (2 * b + 1)) & 1) << b
    return x, y


def _sample(f: Path, d: dict) -> list[dict]:
    """A tile's galaxy sample: from the samples stage when it has run."""
    sp = f.parent.parent / "samples" / f.name
    if sp.exists():
        return json.loads(sp.read_text())["sample"]
    return d.get("sample", [])


def observation_dates(tiles_csv: Path) -> dict[int, tuple[str, str, int]]:
    """Tile -> (first, last, mosaics): a tile can have several VIS mosaics in
    the archive's list (one per visit in the deep-field areas), and the
    survey file keeps only one of them."""
    import csv

    out: dict[int, tuple[str, str, int]] = {}
    if not tiles_csv.exists():
        return out
    for r in csv.DictReader(tiles_csv.open()):
        t = int(r["tile_index"])
        f, la = r["first_observation_date_time"][:10], r["last_observation_date_time"][:10]
        if t in out:
            f0, l0, n = out[t]
            out[t] = (min(f0, f), max(l0, la), n + 1)
        else:
            out[t] = (f, la, 1)
    return out


def pixel_values(
    files: list[Path], dates: dict[int, tuple[str, str, int]] | None = None
) -> tuple[dict[str, dict[int, float]], list[dict]]:
    """Order-PIX_ORDER HEALPix (nested) pixel -> value, per layer; and the
    per-tile summary rows."""
    hp = HEALPix(nside=2**PIX_ORDER, order="nested")
    vals: dict[str, dict[int, float]] = {L.key: {} for L in LAYERS}
    summary = []
    dates = dates or {}
    for f in files:
        d = json.loads(f.read_text())
        t, g = d["tile"], d["grid"]
        if t["tile_index"] in dates:
            first, last, _ = dates[t["tile_index"]]
            t["first_obs"], t["last_obs"] = first, last
        n_all = sum(c[5] for c in d["cells"])
        area = core_area_deg2(t)
        n_gal = sum(
            c[5] for c in d["cells"] if c[2] == 0 and c[4] is not None and c[4] < GAL_MAG_LIMIT
        )
        n_star = sum(c[5] for c in d["cells"] if c[2] is not None and c[2] >= 1)
        summary.append(
            {
                "id": t["tile_index"],
                "ra": round(t["ra"], 5),
                "dec": round(t["dec"], 5),
                "rl": round(t["ra_lo"], 5),
                "dr": round(t["dra"], 5),
                "dl": round(t["dec_lo"], 5),
                "fov": [round(x, 4) for x in t["fov"]],
                "first": t["first_obs"][:10],
                "last": t["last_obs"][:10],
                "n": n_all,
                "ng": n_gal,
                "g": round(n_gal / area, 1),
                "s": round(n_star / area, 1),
                "area": round(area, 5),
                "sample": _sample(f, d)[:SAMPLE_KEEP],
            }
        )
        if not n_all:
            continue
        tv = tile_values(d)
        radius = 0.5 * math.hypot(t["dra"] * math.cos(math.radians(t["dec"])), t["ddec"]) + 0.02
        ip = hp.cone_search_lonlat(t["ra"] * u.deg, t["dec"] * u.deg, radius * u.deg)
        lon, lat = hp.healpix_to_lonlat(ip)
        ra = lon.to_value(u.deg)
        dec = lat.to_value(u.deg)
        fx = ((ra - t["ra_lo"] + 720) % 360) / t["dra"] * g
        fy = (dec - t["dec_lo"]) / t["ddec"] * g
        ok = (fx >= 0) & (fx < g) & (fy >= 0) & (fy < g)
        cx = fx[ok].astype(int)
        cy = fy[ok].astype(int)
        for key, arr in tv.items():
            v = arr[cy, cx]
            m = ~np.isnan(v)
            vals[key].update(zip(ip[ok][m].tolist(), v[m].tolist()))
    return vals, summary


def write_hips(out: Path, layer: Layer, values: dict[int, float]) -> int:
    """Write one HiPS: PNG tiles at orders 0..TILE_ORDER, RGBA, transparent
    where there is no data. Returns the number of tiles written."""
    w = 2**TILE_WIDTH_LOG2
    cmap = _cmap(layer.cmap)
    ipix = np.fromiter(values.keys(), dtype=np.int64)
    v = np.fromiter(values.values(), dtype=float)
    if layer.log:
        x = (np.log10(np.clip(v, 1e-9, None)) - math.log10(layer.vmin)) / (
            math.log10(layer.vmax) - math.log10(layer.vmin)
        )
    else:
        x = (v - layer.vmin) / (layer.vmax - layer.vmin)
    x = np.clip(x, 0, 1)
    written = 0
    # Finest order: one image per parent tile; coarser orders average 2x2.
    level_imgs: dict[int, np.ndarray] = {}  # tile index -> (w, w) float, NaN empty
    parent = ipix >> (2 * TILE_WIDTH_LOG2)
    sub = ipix & ((1 << (2 * TILE_WIDTH_LOG2)) - 1)
    px, py = _deinterleave(sub, TILE_WIDTH_LOG2)
    for tidx in np.unique(parent):
        m = parent == tidx
        img = np.full((w, w), np.nan)
        # Tile pixel layout, checked against Aladin Lite with a synthetic
        # cross: the nested index's EVEN bits run up the image (rows), its
        # odd bits across (columns); _save_tile then flips rows for PNG.
        img[px[m], py[m]] = x[m]
        level_imgs[int(tidx)] = img
    for order in range(TILE_ORDER, -1, -1):
        for tidx, img in level_imgs.items():
            _save_tile(out, order, tidx, img, cmap)
            written += 1
        if order == 0:
            break
        nxt: dict[int, np.ndarray] = {}
        for tidx, img in level_imgs.items():
            par, q = tidx >> 2, tidx & 3
            canvas = nxt.setdefault(par, np.full((w, w), np.nan))
            small = _half(img)
            qx, qy = q & 1, (q >> 1) & 1
            h = w // 2
            canvas[qy * h : (qy + 1) * h, qx * h : (qx + 1) * h] = small
        level_imgs = nxt
    _properties(out, layer)
    return written


def _half(img: np.ndarray) -> np.ndarray:
    h = img.shape[0] // 2
    b = img.reshape(h, 2, h, 2)
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-empty 2x2 blocks stay NaN
        return np.nanmean(b, axis=(1, 3))


def _save_tile(out: Path, order: int, tidx: int, img: np.ndarray, cmap) -> None:
    rgba = (cmap(np.nan_to_num(img)) * 255).astype(np.uint8)
    rgba[..., 3] = np.where(np.isnan(img), 0, 255)
    # Image row 0 is the TOP of a PNG; the HiPS y axis runs upward.
    rgba = rgba[::-1]
    d = out / f"Norder{order}" / f"Dir{(tidx // 10000) * 10000}"
    d.mkdir(parents=True, exist_ok=True)
    Image.fromarray(rgba, "RGBA").save(d / f"Npix{tidx}.png", optimize=True)


def _properties(out: Path, layer: Layer) -> None:
    props = {
        "creator_did": f"ivo://fwh-euclid/P/Euclid/DR1/{layer.key}",
        "obs_title": f"Euclid DR1 wide: {layer.title}",
        "dataproduct_type": "image",
        "hips_version": "1.4",
        "hips_frame": "equatorial",
        "hips_order": str(TILE_ORDER),
        "hips_order_min": "0",
        "hips_tile_width": str(2**TILE_WIDTH_LOG2),
        "hips_tile_format": "png",
        "hips_pixel_scale": f"{58.6323 / 2**PIX_ORDER:.6f}",
        "obs_copyright": "ESA/Euclid/Euclid Consortium, CC BY-NC 3.0 IGO",
    }
    (out / "properties").write_text("".join(f"{k:<20} = {v}\n" for k, v in props.items()))


def _legend(L: Layer) -> dict:
    cm = _cmap(L.cmap)
    stops = [
        "#{:02x}{:02x}{:02x}".format(*tuple(int(255 * c) for c in cm(k / 15)[:3]))
        for k in range(16)
    ]
    if L.log:
        lo, hi = math.log10(L.vmin), math.log10(L.vmax)
        vals = [10 ** (lo + (hi - lo) * k / 4) for k in range(5)]
        ticks = [_short(v) for v in vals]
    elif L.key == "observed":
        from datetime import timedelta

        ticks = [
            (SURVEY_START + timedelta(days=L.vmin + (L.vmax - L.vmin) * k / 4)).strftime("%b %Y")
            for k in range(5)
        ]
    else:
        ticks = [_short(L.vmin + (L.vmax - L.vmin) * k / 4) for k in range(5)]
    return {**L.__dict__, "stops": stops, "ticks": ticks}


def _short(v: float) -> str:
    if v >= 1e6:
        return f"{v / 1e6:.1f}M"
    if v >= 1e3:
        return f"{v / 1e3:.0f}k"
    return f"{v:.0f}"


def _mag(flux_ujy: float | None) -> float | None:
    if not flux_ujy or flux_ujy <= 0:
        return None
    return round(23.9 - 2.5 * math.log10(flux_ujy), 1)


def galaxy_boxes(summary: list[dict]) -> dict[str, list]:
    """The brightest galaxies of every tile, in 10 x 10 degree boxes the
    viewer loads as the view needs them. VIS pixels are 0.1 arcsec."""
    boxes: dict[str, list] = {}
    for t in summary:
        for g in t["sample"]:
            k = f"{int(g['ra'] // 10 * 10) % 360}_{int((g['dec'] + 90) // 10 * 10)}"
            a = g.get("a_px")
            boxes.setdefault(k, []).append(
                [
                    round(g["ra"], 5),
                    round(g["dec"], 5),
                    _mag(g.get("flux")),
                    None if g.get("z") is None else round(g["z"], 2),
                    None if a is None else round(a * 0.1, 1),
                    None if g.get("e") is None else round(g["e"], 2),
                    None if g.get("pa") is None else round(g["pa"]),
                ]
            )
    return boxes


def build(
    survey: Path, out: Path, template: Path | None = None, quasars: Path | None = None
) -> dict:
    files = sorted((survey / "tiles").glob("*.json"))
    vals, summary = pixel_values(files, observation_dates(survey / "tiles.csv"))
    out.mkdir(parents=True, exist_ok=True)
    stats = {}
    for L in LAYERS:
        n = write_hips(out / "hips" / L.key, L, vals[L.key])
        stats[L.key] = {"pixels": len(vals[L.key]), "tiles": n}
    boxes = galaxy_boxes(summary)
    (out / "gal").mkdir(exist_ok=True)
    for k, rows in boxes.items():
        (out / "gal" / f"{k}.json").write_text(json.dumps(rows, separators=(",", ":")))
    slim = [{k: v for k, v in t.items() if k != "sample"} for t in summary]
    (out / "tiles.json").write_text(json.dumps(slim, separators=(",", ":")))
    with_data = [t for t in summary if t["n"]]
    totals = {
        "objects": sum(t["n"] for t in summary),
        "galaxies": sum(t["ng"] for t in summary),
        "area_deg2": sum(t["area"] for t in with_data),
        "tiles": len(summary),
        "tiles_with_data": len(with_data),
        "first": min(t["first"] for t in with_data) if with_data else "",
        "last": max(t["last"] for t in with_data) if with_data else "",
        "galaxy_markers": sum(len(v) for v in boxes.values()),
    }
    site = {
        "totals": totals,
        "gal_mag_limit": GAL_MAG_LIMIT,
        "layers": [_legend(L) for L in LAYERS],
        "gal_boxes": sorted(boxes),
        "jumps": JUMPS,
        "prose_html": "",
    }
    from .prose import prose_html

    site["prose_html"] = prose_html(totals, summary, GAL_MAG_LIMIT)
    if quasars and quasars.exists():
        from .prose import quasar_html

        site["prose_html"] += quasar_html(json.loads(quasars.read_text()))
    tpl = (template or Path(__file__).with_name("site_template.html")).read_text()
    (out / "index.html").write_text(
        tpl.replace("/*__SUMMARY__*/null", json.dumps(site, separators=(",", ":")))
    )
    return {"totals": totals, "layers": stats}


#: Places worth a click. The three Euclid Deep Fields (ESA's published
#: centres) sit inside the wide footprint; the wide views frame each half.
JUMPS = [
    {"name": "Whole north", "ra": 205.0, "dec": 66.0, "fov": 75},
    {"name": "Whole south", "ra": 55.0, "dec": -48.0, "fov": 80},
    {"name": "Deep Field North", "ra": 269.733, "dec": 66.018, "fov": 3},
    {"name": "Deep Field South", "ra": 61.241, "dec": -48.423, "fov": 3},
    {"name": "Deep Field Fornax", "ra": 52.932, "dec": -28.088, "fov": 3},
    {"name": "Large Magellanic Cloud edge", "ra": 80.9, "dec": -62.0, "fov": 20},
    {"name": "All sky", "ra": 180.0, "dec": 0.0, "fov": 200},
]


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(prog="euclidsky.build", description=__doc__.split("\n")[0])
    p.add_argument("--survey", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--quasars", type=Path, help="score JSON from euclidsky.quasar_report")
    a = p.parse_args(argv)
    print(json.dumps(build(a.survey, a.out, quasars=a.quasars), indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
