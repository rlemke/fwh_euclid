# fwh_euclid: ESA Euclid DR1 on the sky

The first Euclid data release (DR1, 12 November 2026) holds about 900 million
objects and 1.7 PB of data. This domain does two things with it, without
downloading it in bulk:

1. **A sky map.** The ESA Euclid Science Archive is asked, one tile at a time,
   to count its own wide-survey catalogue on a 10 x 10 grid laid over each
   tile's core: about 3 arcminutes a cell, ~9,500 tiles. The counts become HiPS
   layers (the IVOA tiled all-sky format) that
   [Aladin Lite](https://aladin.cds.unistra.fr/AladinLite/) draws on the sky as
   seen from Earth, over real imagery. Zoomed in, each tile's brightest galaxies
   appear as ellipses with Euclid's measured size, shape and orientation.
2. **Known quasars, through Euclid's labels.** Every quasar a spectrum has
   already confirmed (Milliquas types Q, A, B, K, N) that falls in a DR1 tile is
   matched to its Euclid object (1 arcsec). The check is what each Euclid label
   calls it: the photometric classifier, the QSO template fit, and Euclid's own
   slitless spectrum, and whether each redshift is right.

## FFL at a glance

```ffl
namespace demo {
    workflow Sky(work_dir: String, site_dir: String) => (index: String) andThen {
        run = euclid.workflows.SkyMap(work_dir = $.work_dir, site_dir = $.site_dir, width = 2)
        yield Sky(index = run.index)
    }
}
```

`SkyMap` runs `euclid.Plan`, then `CountChunk` and `SampleChunk` in fan-outs
capped at `width` chunks, then `BuildSky`. The archive is a shared public
service, so the cap matters more than the speed.

## Command line

```bash
python -m euclidsky.survey  --out DATA/dr1 --workers 8                   # counts
python -m euclidsky.survey  --out DATA/dr1 --workers 8 --stage samples   # galaxy markers
python -m euclidsky.build   --survey DATA/dr1 --out SITE                 # layers + page
python -m euclidsky.quasars --refs refs.csv --out DATA/qso               # cross-match
python -m euclidsky.quasars --refs refs.csv --out DATA/qso --stage reverse
python -m euclidsky.quasar_report --qso DATA/qso --survey DATA/dr1
```

Every stage writes one file per tile and skips tiles already done, so an
interrupted run resumes.

## Things that are not what they seem

- **NaN sorts above every number** in the archive's database. `ORDER BY flux
  DESC` returns NaN fluxes first unless the flux is bounded.
- **The point-like probability is blank** for saturated and NIR-only objects,
  which includes the biggest galaxies. The galaxy markers use the photometric
  class instead, plus a surface-brightness test that rejects saturated stars.
- **The archive's column description of `phz_classification` is reversed.**
  The data-product document (1 = star, 2 = galaxy, 4 = QSO) matches the data.
- **The deep-field areas are in a different table**
  (`mer_catalogue_wide_mode`), with no photometric labels in the `_wide`
  tables. The survey counts them; the quasar check scores only the main survey.
- **Tiles overlap.** Only each tile's core is counted, with the RA width taken
  from the ring's most common neighbour gap.
- **Aladin Lite drew nothing** for 128-pixel HiPS tiles. The tiles are 512
  pixels; the pixel layout inside a tile was fixed by drawing a synthetic cross
  at a known position and keeping the one orientation of eight that put it
  there (`test_hips_pixel_lands_where_its_sky_position_is`).
- **The archive was flaky on release day.** It sometimes answered "Unknown
  table" for tables that exist. That error is retried as transient.

Data: ESA/Euclid/Euclid Consortium, CC BY-NC 3.0 IGO. Reference quasars:
Milliquas (Flesch), via NASA HEASARC.
