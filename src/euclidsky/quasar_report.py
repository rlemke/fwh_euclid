"""Score the quasar cross-match: what Euclid DR1 calls quasars that are known.

Only references inside the wide-survey catalogue are scored, and only where
Euclid actually has data at that spot: a tile can have holes (masks around
bright stars, gaps between exposures), and a quasar in a hole is a fact
about coverage, not about classification. Coverage is read from the survey's
own counts: the reference's ~3 arcmin cell must hold at least a fifth of the
tile's median cell count.
"""

from __future__ import annotations

import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path

from ._tap import num

SURVEY_CAT = "dr1.mer_catalogue_wide_survey"
MAG_BINS = [
    (0, 19, "brighter than 19"),
    (19, 20, "19–20"),
    (20, 21, "20–21"),
    (21, 22, "21–22"),
    (22, 99, "fainter than 22"),
]


def _covered(survey: Path, tile: int, ra: float, dec: float, cache: dict) -> bool | None:
    if tile not in cache:
        p = survey / "tiles" / f"{tile}.json"
        if not p.exists():
            cache[tile] = None
        else:
            d = json.loads(p.read_text())
            per = Counter()
            for c in d["cells"]:
                per[(c[0], c[1])] += c[5]
            med = statistics.median(per.values()) if per else 0
            cache[tile] = (d["tile"], d["grid"], per, med)
    v = cache[tile]
    if v is None:
        return None
    t, g, per, med = v
    cx = int(((ra - t["ra_lo"] + 720) % 360) / t["dra"] * g)
    cy = int((dec - t["dec_lo"]) / t["ddec"] * g)
    return per.get((cx, cy), 0) >= 0.2 * med


def load(qso: Path, survey: Path) -> list[dict]:
    lab: dict[str, dict] = {}
    for p in (qso / "labels").glob("*.json"):
        lab.update(json.loads(p.read_text()))
    rows = []
    cache: dict = {}
    for p in sorted((qso / "xmatch").glob("*.json")):
        for o in json.loads(p.read_text()):
            if o["catalogue"] != SURVEY_CAT:
                continue
            r = o["ref"]
            ra, dec = float(r["ra"]), float(r["dec"])
            cov = _covered(survey, o["tile"], ra, dec, cache)
            m = o["match"]
            L = lab.get(m["object_id"], {}) if m else {}
            rows.append(
                {
                    "type": r["type"],
                    "z": num(r["z"]),
                    "rmag": num(r["rmag"]),
                    "covered": cov,
                    "matched": bool(m),
                    "labelled": bool(L),
                    "phz_class": None
                    if L.get("phz_classification") in (None, "")
                    else int(float(L["phz_classification"])),
                    "qso_prob": num(L.get("phz_qso_prob")),
                    "qso_threshold": num(L.get("qso_threshold")),
                    "qso_fit": "qso_fit_z" in L,
                    "qso_fit_z": num(L.get("qso_fit_z")),
                    "phz_median": num(L.get("phz_median")),
                    "spe_class": L.get("spe_class") or None,
                    "spe_qso_z": num(L.get("spe_qso_z")),
                    "spe_gal_z": num(L.get("spe_gal_z")),
                }
            )
    return rows


def _frac(a: int, b: int) -> float | None:
    return round(a / b, 4) if b else None


def _dz(a: float | None, b: float | None) -> float | None:
    if a is None or b is None:
        return None
    return abs(a - b) / (1 + b)


def score(rows: list[dict]) -> dict:
    q = [r for r in rows if r["type"].startswith("Q") and r["covered"]]
    det = [r for r in q if r["matched"] and r["labelled"]]
    spec = [r for r in det if r["spe_class"]]
    out = {
        "references_in_survey": sum(1 for r in rows if r["type"].startswith("Q")),
        "covered": len(q),
        "detected": sum(1 for r in q if r["matched"]),
        "labelled": len(det),
        # Negative classes mean "not classified", and -1 has every bit set.
        "phz_qso_flag": sum(1 for r in det if (r["phz_class"] or 0) > 0 and r["phz_class"] & 4),
        "phz_class": dict(
            Counter("none" if r["phz_class"] is None else str(r["phz_class"]) for r in det)
        ),
        "qso_prob_positive": sum(1 for r in det if (r["qso_prob"] or 0) > 0),
        "qso_thresholds": dict(Counter(str(r["qso_threshold"]) for r in det)),
        "qso_fit": sum(1 for r in det if r["qso_fit"]),
        "with_spectrum": len(spec),
        "spe_class": dict(Counter(r["spe_class"] for r in spec)),
    }

    # Redshifts. A quasar's own redshift from Euclid's spectrum is right when
    # within 1% in (1+z); from photometry, within 10%.
    def acc(vals: list[float | None], tol: float) -> dict:
        v = [x for x in vals if x is not None]
        return {
            "n": len(v),
            "within": sum(1 for x in v if x < tol),
            "frac": _frac(sum(1 for x in v if x < tol), len(v)),
        }

    zr = [r for r in det if r["z"] is not None]
    out["z_spec_qso_template"] = acc(
        [_dz(r["spe_qso_z"], r["z"]) for r in zr if r["spe_class"]], 0.01
    )
    out["z_spec_reported"] = acc(
        [
            _dz(r["spe_qso_z"] if r["spe_class"] == "qso" else r["spe_gal_z"], r["z"])
            for r in zr
            if r["spe_class"]
        ],
        0.01,
    )
    out["z_spec_reported_by_class"] = {
        c: acc(
            [
                _dz(r["spe_qso_z"] if c == "qso" else r["spe_gal_z"], r["z"])
                for r in zr
                if r["spe_class"] == c
            ],
            0.01,
        )
        for c in ("qso", "galaxy", "star")
    }
    out["z_photo_qso_fit"] = acc([_dz(r["qso_fit_z"], r["z"]) for r in zr if r["qso_fit"]], 0.1)
    out["z_photo_galaxy"] = acc([_dz(r["phz_median"], r["z"]) for r in zr], 0.1)
    by_mag = []
    for lo, hi, name in MAG_BINS:
        sel = [r for r in q if r["rmag"] is not None and lo <= r["rmag"] < hi]
        d = [r for r in sel if r["matched"] and r["labelled"]]
        s = [r for r in d if r["spe_class"]]
        by_mag.append(
            {
                "bin": name,
                "known": len(sel),
                "detected": _frac(sum(1 for r in sel if r["matched"]), len(sel)),
                "spectrum": _frac(len(s), len(d)),
                "spec_says_qso": _frac(sum(1 for r in s if r["spe_class"] == "qso"), len(s)),
                "qso_template_z_ok": _frac(
                    sum(1 for r in s if (_dz(r["spe_qso_z"], r["z"]) or 9) < 0.01), len(s)
                ),
            }
        )
    out["by_mag"] = by_mag
    by_z = defaultdict(lambda: [0, 0])
    for r in det:
        if r["spe_class"] and r["z"] is not None:
            k = min(int(r["z"] // 0.5) * 0.5, 3.0)
            by_z[k][0] += 1
            by_z[k][1] += r["spe_class"] == "qso"
    out["spec_says_qso_by_z"] = [
        {"z": f"{k:.1f}–{k + 0.5:.1f}" if k < 3 else "3.0+", "n": n, "frac": _frac(s, n)}
        for k, (n, s) in sorted(by_z.items())
    ]
    return out


def reverse(qso: Path, refs_csv: Path) -> dict:
    """Of the objects Euclid's spectra call quasars, in the tiles where the
    reference is densest: how many are known (any confirmed active nucleus
    in the reference within 1")?"""
    import csv
    import math

    ref = defaultdict(list)
    for r in csv.DictReader(refs_csv.open()):
        ref[int(r["tile"])].append((float(r["ra"]), float(r["dec"]), r["type"]))
    n = known = cand = 0
    bright = bright_starlike = faint = faint_known = 0
    for p in (qso / "reverse").glob("*.json"):
        d = json.loads(p.read_text())
        refs = ref.get(d["tile"], [])
        for o in d["rows"]:
            ra, dec = float(o["right_ascension"]), float(o["declination"])
            n += 1
            best = None
            for rra, rdec, typ in refs:
                dd = math.hypot((rra - ra) * math.cos(math.radians(dec)), rdec - dec) * 3600
                if dd < 1.0:
                    best = typ
                    break
            is_known = bool(best and best[:1] in "QABKN")
            if is_known:
                known += 1
            elif best:
                cand += 1
            f = num(o.get("flux_detection_total"))
            mag = -2.5 * math.log10(f) + 23.9 if f and f > 0 else None
            pl = num(o.get("point_like_prob"))
            if mag is not None and mag < 20:
                bright += 1
                bright_starlike += pl is None or pl >= 0.5
            elif mag is not None and mag >= 22:
                faint += 1
                faint_known += is_known
    return {
        "tiles": len(list((qso / "reverse").glob("*.json"))),
        "spectral_qsos": n,
        "known": known,
        "reference_candidates_only": cand,
        "known_frac": _frac(known, n),
        "bright": bright,
        "bright_starlike": bright_starlike,
        "faint": faint,
        "faint_known": faint_known,
    }


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser(prog="euclidsky.quasar_report")
    ap.add_argument("--qso", type=Path, required=True)
    ap.add_argument("--survey", type=Path, required=True)
    ap.add_argument("--out", type=Path)
    a = ap.parse_args()
    res = score(load(a.qso, a.survey))
    refs = a.qso / "mq_in_tiles.csv"
    if (a.qso / "reverse").exists() and refs.exists():
        res["reverse"] = reverse(a.qso, refs)
    text = json.dumps(res, indent=1)
    if a.out:
        a.out.write_text(text)
    print(text)


if __name__ == "__main__":
    main()
