# 02_GetRubinVisits — LSSTCam Visit Tables with Tract/Patch Indices

**Author:** Sylvie Dagoret-Campagne — IJCLab / IN2P3 / CNRS — Université Paris-Saclay  
**Created:** 2026-03-25  
**Last updated:** 2026-06-24  
**Context:** Rubin/LSST DP2 — Building visit tables enriched with Butler skymap
tract/patch indices for use in downstream photometric analyses.

---

## Overview

This directory contains two complementary notebooks that both produce the same
kind of output — a **visit table enriched with `(tract, patch)` indices** — but
using two different data sources:

| Source | Notebook | Output directory |
|--------|----------|-----------------|
| Butler registry (`main` / `LSSTCam/defaults`) | `2026-03-25_FindObservationsInButlerRegistryInTractPatch.ipynb` | `data_DP2_VisitsAndTracts/` |
| ConsDB (`cdb_lsstcam.visit1`) | `2026-03-26_DP2_ConstDB_Butler_LSSTCam_VisitsTractPatch.ipynb` | `data_DP2_ConstDB_VisitsAndTracts/` |

The enriched visit tables are the primary input of notebook series
`01_askSimbad_stars` (NB02: cross-match with stable stars) and
`03_VisualiseCoadds` (patch selection for deepCoadd loading).

---

## DDF fields covered

| Key | RA (deg) | Dec (deg) |
|-----|----------|-----------|
| COSMOS    | 150.119 | +2.206  |
| ECDFS     |  53.125 | −28.100 |
| ELAIS-S1  |   9.450 | −44.000 |
| XMM-LSS   |  35.708 | −4.750  |
| EDFS-a    |  58.900 | −49.315 |
| EDFS-b    |  63.600 | −47.600 |
| EDFS      |  61.240 | −48.423 |
| M49       | 187.400 | +8.000  |

Search / cone radius: **1.8 deg** from each field centre.

---

## Notebooks

### `2026-03-25_FindObservationsInButlerRegistryInTractPatch.ipynb`

**Goal:** Query the Butler registry directly for all LSSTCam science exposures
and associate each visit pointing with its skymap `(tract, patch)`.

**Butler configuration:**

| Parameter | Value |
|-----------|-------|
| `repo` | `main` |
| `collection` | `LSSTCam/defaults` |
| `skymap` | `lsst_cells_v2` |
| `date_start` | 20250415 |
| `date_stop` | 20260630 |

**Key steps:**

1. Query `registry.queryDimensionRecords("exposure", where=WHERE_CLAUSE_DATE)`
   for all LSSTCam science exposures in the date range.
2. Convert Butler `Timespan` JD values to MJD and ISO-T strings via
   `astropy.time.Time`.
3. Filter to `observation_type == "science"` only.
4. For each visit pointing `(ra, dec)`, call `skymap.findTract()` then
   `tract_info.findPatch()` to obtain `tract`, `patch` (sequential index),
   `patch_str` (2D index string), and the tract bounding box corners in both
   pixel and sky coordinates.
5. Save the result as a Parquet file named
   `visitTable-{id_min}-{id_max}_N{nvisits}_WithTracts.parquet`.

**Key function:** `getLostOfVisits(registry, where_clause_date)` —
iterates over dimension records, builds rows with timing, pointing, filter,
and zenith angle fields, and returns a filtered `DataFrame`.

**Output columns (selected):**

| Column | Description |
|--------|-------------|
| `id` | Butler visit/exposure ID |
| `mjd` | MJD of exposure start |
| `filter` | Physical filter name (e.g. `LSSTCam-r`) |
| `band` | Derived band letter (e.g. `r`) |
| `ra`, `dec` | Tracking pointing (deg) |
| `zenith_angle` | Zenith angle (deg) |
| `tract`, `patch` | Skymap indices (`Int64`) |
| `patch_str` | 2D patch index string (e.g. `"5,3"`) |
| `tract_bbox` | Pixel bounding box tuple `(xmin, ymin, xmax, ymax)` |
| `tract_ra_corners`, `tract_dec_corners` | Corner RA/Dec lists |

---

### `2026-03-26_DP2_ConstDB_Butler_LSSTCam_VisitsTractPatch.ipynb`

**Goal:** Same tract/patch enrichment, but using **ConsDB** as the visit
source instead of the Butler registry. This provides access to additional
ConsDB-only columns (e.g. `exp_time`, `s_ra`, `s_dec`, `physical_filter`
from `cdb_lsstcam.visit1`).

**ConsDB configuration:**

| Parameter | Value |
|-----------|-------|
| URL | `http://consdb-pq.consdb:8080/consdb` |
| Table | `cdb_lsstcam.visit1` |
| Filter date | `day_obs >= 20250415` |

**Key steps:**

1. Connect to ConsDB via `lsst.summit.utils.ConsDbClient`.
2. Query `cdb_lsstcam.visit1` for all visits since 2025-04-15 and
   inner-join with itself (following the parent notebook pattern).
3. **Clean visits:**
   - Remove engineering/pinhole filters (`other`, `none`, `other:pinhole`, `ph_5`).
   - Drop rows with missing pointing coordinates (`s_ra`, `s_dec`).
   - Select only science exposures with `20 s ≤ exp_time ≤ 40 s`.
4. Initialise Butler (`main` / `LSSTCam/defaults`) and load the skymap.
5. Apply the same `get_tract_patch()` function row-by-row (with `raname="s_ra"`,
   `decname="s_dec"`) to enrich with `tract`, `patch`, `patch_str`,
   `tract_bbox`, and sky corner columns.
6. Save as
   `constDbVisitTable-{id_min}-{id_max}_N{nvisits}_WithTracts.parquet`.

---

## Shared utility: `get_tract_patch()`

Both notebooks use the same helper function:

```python
def get_tract_patch(row, skymap, raname="ra", decname="dec"):
    sp         = SpherePoint(row[raname], row[decname], degrees)
    tract_info = skymap.findTract(sp)
    patch_info = tract_info.findPatch(sp)
    wcs        = tract_info.getWcs()
    # → returns tract, patch (sequential), patch_str (2D),
    #   tract_bbox (pixel), tract_ra/dec_corners (sky)
```

This converts a pointing `(RA, Dec)` to a `lsst_cells_v2` tract/patch pair,
also computing the tract sky footprint corners for downstream spatial selection.

---

## Output files

| File pattern | Source notebook | Description |
|-------------|-----------------|-------------|
| `data_DP2_VisitsAndTracts/visitTable-{min}-{max}_N{N}_NoTracts.csv` | NB 2026-03-25 | Raw visit list before tract enrichment |
| `data_DP2_VisitsAndTracts/visitTable-{min}-{max}_N{N}_WithTracts.parquet` | NB 2026-03-25 | Butler registry visits + tract/patch |
| `data_DP2_ConstDB_VisitsAndTracts/constDbVisitTable-{min}-{max}_N{N}_WithTracts.parquet` | NB 2026-03-26 | ConsDB visits + tract/patch |

The `_WithTracts.parquet` files are the canonical outputs consumed by
downstream notebooks.

---

## Software environment

| Item | Value |
|------|-------|
| Python | 3.12 (LSST kernel) |
| Key LSST packages | `lsst.daf.butler`, `lsst.geom`, `lsst.skymap`, `lsst.summit.utils` |
| Key Python packages | `astropy`, `pandas`, `numpy`, `matplotlib` |
| Execution environment | USDF / RSP JupyterLab |
