"""The page's explanatory text, with its numbers filled in from the build."""

from __future__ import annotations

import statistics


def _f(x: float) -> str:
    return f"{x:,.0f}"


def prose_html(totals: dict, summary: list[dict], mag: int) -> str:
    rows = [t for t in summary if t["n"]]
    empty = len(summary) - len(rows)
    g = sorted(t["g"] for t in rows)
    med = statistics.median(g) if g else 0
    p10 = g[len(g) // 10] if g else 0
    p90 = g[(len(g) * 9) // 10] if g else 0
    s = sorted(t["s"] for t in rows)
    smed = statistics.median(s) if s else 0
    smax = s[-1] if s else 0
    per_arcmin = med / 3600
    return f"""
<h2>What you are looking at</h2>
<p>The map is the inside of the celestial sphere, as seen from Earth: north is up, east is to the
left, and the background is a photograph of the real night sky. Euclid looks at patches of it far
from the Milky Way's dusty plane, where the view out of our galaxy is clearest. Those patches are
the coloured areas. DR1 covers about {_f(totals["area_deg2"])} square degrees in
{_f(totals["tiles_with_data"])} tiles of roughly half a degree, about the size of the full Moon.</p>
<p>The default colour is the number of galaxies brighter than magnitude {mag} per square degree.
A typical tile has {_f(med)} of them per square degree (the middle 80% of tiles range from
{_f(p10)} to {_f(p90)}): about {per_arcmin:.0f} in a patch of sky one arcminute across, the size of
a grain of rice held at arm's length. The colour is smooth because the universe is nearly uniform
on these scales. The variations that remain come mostly from the survey, not from cosmic
structure: shallower exposures, bright stars masking their surroundings, and dust.</p>
<p>Zoom in below about 8 degrees and the brightest galaxies of each tile appear as ellipses with
Euclid's measured size and shape, coloured by their photometric redshift (blue nearest, then green,
yellow, red farther than about z = 1.2). Click one for its numbers. These are a sample, 30 per tile
({_f(totals["galaxy_markers"])} in all), not the {totals["galaxies"] / 1e6:,.0f} million galaxies
counted.</p>

<h2>How it was made</h2>
<ul>
<li>Nothing was downloaded in bulk. The ESA Euclid Science Archive was asked, one tile at a time,
to count its own catalogue (<code>dr1.mer_catalogue_wide_survey</code>, {_f(totals["objects"])}
objects after removing those flagged spurious) on a 10 x 10 grid laid over each tile, about 3
arcminutes a cell.</li>
<li>Tiles overlap their neighbours a little. Only each tile's core was counted, so every object is
counted once.</li>
<li>A galaxy here is an object Euclid measured as extended (point-like probability below 0.5) and
brighter than magnitude {mag} in its detection image. The magnitude cut keeps the count
comparable between deep and shallow tiles. Point-like objects are mostly Milky Way stars: their
density is {_f(smed)} per square degree in a typical tile and up to {_f(smax)} near the plane of
the Galaxy and the Magellanic Clouds.</li>
<li>{empty} tiles in the archive's tile list have no objects in the wide catalogue. They are shown
in their outline only.</li>
</ul>

<h2>Caveats</h2>
<ul>
<li>This is DR1-Foundation, released on 12 November 2026. Some catalogue fields are
documented inconsistently in the archive: the archive's column description of the star/galaxy
class flag gives star and galaxy the opposite way round from the data product document. The
document's coding matches the data and is the one used here.</li>
<li>The galaxy layer does not correct for depth, masking or dust, so it is a map of what Euclid
detected, not of how galaxies are distributed in space.</li>
</ul>
"""


def _pct(x: float | None) -> str:
    return "–" if x is None else f"{100 * x:.0f}%"


def _bars(items: list[tuple[str, int, str]], total: int) -> str:
    """A horizontal bar chart (one series, so no legend), share of ``total``."""
    rows = []
    for name, n, note in items:
        w = 100 * n / total if total else 0
        rows.append(
            f'<div class="qbar" title="{name}: {n:,} of {total:,} ({w:.1f}%)">'
            f'<span class="qlab">{name}</span>'
            f'<span class="qtrack"><span class="qfill" style="width:{w:.1f}%"></span></span>'
            f'<span class="qval">{w:.0f}%</span><span class="qnote">{note}</span></div>'
        )
    return '<div class="qbars">' + "".join(rows) + "</div>"


def quasar_html(s: dict) -> str:
    """The quasar section, from ``quasar_report.score``."""
    det = s["labelled"]
    spec = s["with_spectrum"]
    sc = s["spe_class"]
    zq = s["z_spec_qso_template"]
    zr = s["z_spec_reported"]
    zg = s["z_spec_reported_by_class"].get("galaxy", {})
    fit = s["z_photo_qso_fit"]
    rows = "".join(
        f"<tr><td>{b['bin']}</td><td>{b['known']:,}</td><td>{_pct(b['detected'])}</td>"
        f"<td>{_pct(b['spectrum'])}</td><td>{_pct(b['spec_says_qso'])}</td>"
        f"<td>{_pct(b['qso_template_z_ok'])}</td></tr>"
        for b in s["by_mag"]
        if b["known"]
    )
    zrows = "".join(
        f"<tr><td>{b['z']}</td><td>{b['n']:,}</td><td>{_pct(b['frac'])}</td></tr>"
        for b in s["spec_says_qso_by_z"]
    )
    rev = s.get("reverse")
    rev_html = ""
    if rev and rev["spectral_qsos"]:
        rev_html = f"""
<p><b>The other direction.</b> In the {rev["tiles"]} tiles where the reference is densest,
Euclid's spectra call {rev["spectral_qsos"]:,} objects quasars. Only
{rev["known"]:,} ({_pct(rev["known_frac"])}) are confirmed active nuclei in the reference.
That is a lower bound on how many are real, because no reference has taken a spectrum
of every faint object. It does show the label cannot be read as "a quasar" without
checking.</p>"""
    bars = _bars(
        [
            ("galaxy", sc.get("galaxy", 0), "the redshift reported is the galaxy template's"),
            ("quasar", sc.get("qso", 0), ""),
            ("star", sc.get("star", 0), ""),
        ],
        spec,
    )
    return f"""
<h2>Does Euclid recognise quasars that are already known?</h2>
<p>ESA's release counts "over 60 million quasar candidates". No table in the public archive
marks such a set, so this checks the labels the archive does carry against
{s["covered"]:,} quasars whose spectra other telescopes have already taken (the Million
Quasars catalogue, type-1 quasars only), at places where Euclid has data.
Euclid detects {_pct(s["detected"] / s["covered"] if s["covered"] else None)} of them
within 1 arcsecond.</p>
<ul>
<li><b>The photometric quasar classifier is not switched on in DR1.</b> All {det:,} detected
quasars have a quasar probability of 0 against a threshold of 1.0, and the quasar bit of
the class flag is set for {s["phz_qso_flag"]:,} of them. The classifier calls
{_pct(s["phz_class"].get("2", 0) / det if det else None)} of them galaxies.</li>
<li><b>A quasar template is fitted to the photometry</b> of {_pct(s["qso_fit"] / det if det else None)}
of them. Its redshift is within 10% (in 1+z) for {_pct(fit["frac"])} of {fit["n"]:,}.</li>
<li><b>Euclid's own spectrum</b> exists for {spec:,} of them. Its classifier calls them:</li>
</ul>
{bars}
<p>Fitting the quasar template to those same spectra gives the right redshift (within 1%
in 1+z) for {_pct(zq["frac"])} of them. The redshift Euclid actually reports, from whichever
template its classifier picked, is right for {_pct(zr["frac"])}. For the quasars it calls
galaxies, the reported redshift is right for {_pct(zg.get("frac"))}. The spectra see the
quasar; the classifier picks the wrong template.</p>
{rev_html}
<table>
<thead><tr><th>Known quasar, r magnitude</th><th>Known</th><th>Detected</th><th>Has a Euclid spectrum</th>
<th>Spectrum says quasar</th><th>Quasar-template z right</th></tr></thead>
<tbody>{rows}</tbody></table>
<p style="margin-top:12px"></p>
<table>
<thead><tr><th>Known redshift</th><th>With a Euclid spectrum</th><th>Spectrum says quasar</th></tr></thead>
<tbody>{zrows}</tbody></table>
<p class="hint">Euclid's slitless spectra cover 1.2–1.9 µm. A quasar's broad lines fall in
that window only at some redshifts, which is why the classification changes with z.
Reference: Milliquas (Flesch), entries confirmed by a spectrum, via NASA HEASARC.
Matching: the nearest Euclid object within 1 arcsec.</p>
"""
