# FFL examples — `euclid`

Complete, compile-checked scenarios against this domain's own facets. Every
fenced `ffl` block below compiles:

```bash
fw ffl compile --check --primary <block>.ffl \
  --library ~/fw_handlers/fwh_euclid/src/euclidsky/ffl/euclid.ffl \
  --library facetwork/ffl/provenance.ffl
```

Example 6 also needs the shared publisher's FFL from fwh_census_us
(`census_publish.ffl` and `census_operations.ffl`).

The facets in play:

| Namespace | Facet | Returns |
|---|---|---|
| `euclid` | `Plan(work_dir, chunk_size)` | `chunks, tiles` |
| `euclid` | `CountChunk(work_dir, tiles, workers)` | `counted, failed` |
| `euclid` | `SampleChunk(work_dir, tiles, workers)` | `sampled, failed` |
| `euclid` | `BuildSky(work_dir, site_dir)` | `index, objects, galaxies, area_deg2` |

And the workflows built from them: `euclid.workflows.CountAll` and `SampleAll`
(a capped `foreach` over chunks) and `SkyMap` (plan, count, sample, build).

Every step writes one file per tile under `work_dir` and skips tiles already
done, so a retried or repeated step does only what is missing. ⚠️ The ESA
Euclid Science Archive is a shared public service: keep the fan-out narrow.

---

## 1. The minimal workflow

The whole page in one call. `width` is how many chunks are counted at once.

```ffl
namespace ex.minimal {
    workflow Sky(work_dir: String, site_dir: String) => (index: String, objects: Long) andThen {
        run = euclid.workflows.SkyMap(work_dir = $.work_dir, site_dir = $.site_dir, width = 2)
        yield Sky(index = run.index, objects = run.objects)
    }
}
```

## 2. Density only, no galaxy markers

Skip the sampling stage. The build reads what the counts wrote to `work_dir`,
a dependency no value carries, so it is ordered with `after`.

```ffl
namespace ex.counts {
    workflow DensityOnly(work_dir: String, site_dir: String) => (index: String, area_deg2: Double) andThen {
        plan = euclid.Plan(work_dir = $.work_dir)
        counts = euclid.workflows.CountAll(work_dir = $.work_dir, chunks = plan.chunks, width = 2)
        site = euclid.BuildSky(work_dir = $.work_dir, site_dir = $.site_dir) after counts
        yield DensityOnly(index = site.index, area_deg2 = site.area_deg2)
    }
}
```

## 3. Gentle on the archive

Small chunks counted one at a time (`limit 1` serialises), each with a single
query thread. Slower, and it adds almost nothing to the archive's load.

```ffl
namespace ex.gentle {
    workflow OneAtATime(work_dir: String) => (counted: Json, tiles: Long) andThen {
        plan = euclid.Plan(work_dir = $.work_dir, chunk_size = 50)
        counts = Counts(work_dir = $.work_dir, chunks = plan.chunks)
        yield OneAtATime(counted = counts.counted, tiles = plan.tiles)
    }

    facet Counts(work_dir: String, chunks: Json) => (counted: Json, failed: Json)
        andThen foreach chunk in $.chunks limit 1 {
            c = euclid.CountChunk(work_dir = $.work_dir, tiles = $.chunk, workers = 1)
            yield Counts(counted += c.counted, failed += c.failed)
        }
}
```

## 4. Do not publish an empty page

On release day the archive sometimes answered "Unknown table" for tables that
exist. A run that counted nothing builds a valid but empty page; here that
case returns a status instead of an index.

```ffl
namespace ex.gate {
    workflow PageIfData(work_dir: String, site_dir: String) => (index: String, status: String) andThen {
        run = euclid.workflows.SkyMap(work_dir = $.work_dir, site_dir = $.site_dir) andThen when {
            case $.objects > 0 => {
                yield PageIfData(index = $.index, status = "ok")
            }
            case _ => {
                yield PageIfData(index = "", status = "empty")
            }
        }
    }
}
```

## 5. Record what produced the page

`fw.provenance.Manifest` writes the software versions and a digest of every
output. The tile counts were not fetched through `fw.http`, so they carry no
provenance sidecar; the manifest reports them in `incomplete` rather than
leaving them out.

```ffl
namespace ex.provenance {
    workflow SkyWithManifest(work_dir: String, site_dir: String) => (index: String, complete: Boolean, incomplete: [String]) andThen {
        run = euclid.workflows.SkyMap(work_dir = $.work_dir, site_dir = $.site_dir)
        man = fw.provenance.Manifest(outputs = $.site_dir,
                                     inputs = [$.work_dir ++ "/tiles.csv"],
                                     dest = $.site_dir ++ "/provenance.json") after run
        yield SkyWithManifest(index = run.index, complete = man.complete, incomplete = man.incomplete)
    }
}
```

## 6. Build and publish (cross-domain)

The shared Pages publisher from fwh_census_us replaces only the `dests`
subtree and merges the gallery landing. ⚠️ `labels` sets the landing link
text: pass the existing label on a republish, or the link is renamed.

```ffl
namespace ex.publish {
    use census.Publish

    workflow BuildAndPublish(work_dir: String, site_dir: String) => (pages_url: String, commit: String) andThen {
        run = euclid.workflows.SkyMap(work_dir = $.work_dir, site_dir = $.site_dir)
        pub = PublishWebBundle(repo = "rlemke/facetwork-maps", prefixes = [$.site_dir],
                               dests = ["space/euclid-dr1"],
                               labels = ["Euclid DR1: galaxies on the sky"]) after run
        yield BuildAndPublish(pages_url = pub.pages_url, commit = pub.commit)
    }
}
```

---

## Syntax cheat sheet

| | |
|---|---|
| step | `name = ns.Facet(arg = value)` |
| reference | `step.field`, never a bare step name |
| container attributes | `$.param` (this block), `$$.param` (one up) |
| string concatenation | `$.work_dir ++ "/tiles.csv"` |
| fan-out | `andThen foreach x in $.list limit $.width { … yield F(out += step.field) }` |
| branch | `step = F(…) andThen when { case $.n > 0 => { … } case _ => { … } }` |
| ordering without a value | `b = F(…) after a` |
