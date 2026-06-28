# RubinLSSTStableStarsinDDF — Notebook Series Overview

**Author:** Sylvie Dagoret-Campagne — IJCLab / IN2P3 / CNRS, Université Paris-Saclay  
**Project:** Identification and characterisation of photometrically stable stars
in Rubin/LSST Deep Drilling Fields (DDFs) for use as calibration anchors.  
**Execution environment:** Rubin Science Platform (RSP) at USDF, or local Mac
(`conda_py313`, Python 3.13) depending on the series.

---

## Scientific context

Photometrically stable stars with known spectral types are essential calibration
anchors for the Rubin/LSST survey. In the Deep Drilling Fields (COSMOS,
XMM-LSS, ELAIS-S1, ECDFS, EDFS), a large number of repeated visits per field
makes it possible to characterise atmospheric effects (DCR, PWV, airmass) and
instrument systematics at the level of individual visits. This project builds a
pipeline that goes from an astronomical database query (SIMBAD) all the way to
per-star multi-band photometric light curves extracted from the Rubin Gen3
Butler and displayed on deepCoadd images.

---

## Pipeline overview

The five subdirectories form a **sequential pipeline**. Each series produces
output files consumed by the next one.

```
01_askSimbad_stars
        │  stable-star master catalogue + per-DDF target lists
        ▼
02_GetRubinVisits
        │  visit table (tract/patch enriched)
        ▼
03_VisualiseCoadds ◄── (also uses output of 04_FindObjects)
        │  deepCoadd images + stable-star overlays in Firefly
        │
04_FindObjects
        │  Simbad × LSSTCam object cross-match + coadd-level photometry
        ▼
05_FindSources
           per-visit source light curves (psfFlux vs. MJD)
```

---

## Subdirectories

### `01_askSimbad_stars/`

**Topic:** Query SIMBAD for photometrically stable, non-variable stars in
the DDFs and build a master target catalogue.

**What it does:**
- Cone-searches SIMBAD around each DDF centre (radius 1.5 deg, V ≤ 20)
  for non-variable stellar objects (`otype = *`) with known spectral types.
- Parses MK spectral types (temperature class, subtype, luminosity class)
  and produces diagnostic distributions (sky maps, V-mag histograms,
  spectral-type bar charts, parallax distributions).
- Cross-matches the master catalogue against the Rubin visit table
  (angular separation < 1.75 deg from pointing) to count per-band visits
  per star and compute airmass time series.
- Retrieves PanSTARRS DR1 HiPS cutouts via CDS HiPS2FITS for visual QC.
- Exports per-DDF target lists in CSV format for upload to the RSP Portal.

**Key output:** `data_SIMBAD_04/stable_stars_DDF<field>.csv` — one file per
DDF, consumed by `04_FindObjects` and `05_FindSources`.

**Execution environment:** local Mac (`conda_py313`)  
**README:** [`01_askSimbad_stars/README.md`](01_askSimbad_stars/README.md)

---

### `02_GetRubinVisits/`

**Topic:** Build complete LSSTCam visit tables enriched with Butler skymap
tract/patch indices.

**What it does:**
- Provides two complementary approaches to obtain the visit list:
  (1) querying the Butler registry directly (`main` / `LSSTCam/defaults`),
  (2) querying ConsDB (`cdb_lsstcam.visit1`) via `ConsDbClient`.
- For each visit pointing `(RA, Dec)`, resolves the `lsst_cells_v2` tract
  and patch using `skymap.findTract()` / `findPatch()`, and stores the tract
  bounding box and sky-corner coordinates.
- ConsDB variant applies additional cleaning: removes engineering filters
  and restricts to science exposures with `20 s ≤ exp_time ≤ 40 s`.
- Outputs Parquet files named
  `visitTable-{id_min}-{id_max}_N{N}_WithTracts.parquet` (Butler) and
  `constDbVisitTable-..._WithTracts.parquet` (ConsDB).

**Key output:** `_WithTracts.parquet` visit tables — used by
`01_askSimbad_stars` (NB02 visit cross-match) and `05_FindSources`
(MJD timestamp merge).

**Execution environment:** RSP / USDF (LSST kernel)  
**README:** [`02_GetRubinVisits/README.md`](02_GetRubinVisits/README.md)

---

### `03_VisualiseCoadds/usdf_butler/`

**Topic:** Display LSSTCam deepCoadd images in Firefly with stable-star
overlays annotated with SIMBAD name and spectral type.

**What it does:**
- Provides three notebooks with increasing complexity:
  (1) mosaic of deepCoadd patches displayed in Firefly,
  (2) same mosaic in Matplotlib with WCS projection,
  (3) deepCoadd display with overlay of photometrically stable stars.
- **Patch-coverage analysis:** for each DDF, counts the distinct
  `(tract, patch)` pairs covering the stable stars and applies a
  `npatchmax` decision rule:
  - `n_patches > npatchmax` → load patches in the single selected band
    (`BANDSEL`, default `r`),
  - `n_patches ≤ npatchmax` → load all 6 bands (u g r i z y).
- Converts star `(RA, Dec)` coordinates to image pixel coordinates via
  `wcs.skyToPixel()` and draws markers with `afw_display.dot()` inside
  a `Buffering()` context.
- Each Firefly frame title includes the field name, band, tract, patch
  number, and star count. Markers are annotated with the SIMBAD identifier
  and spectral type.
- Includes a Matplotlib sanity-check preview saved as PDF + PNG.

**Key input:** `data_DEEPCOADDSOBJ_03_in/targets_matched_lsst_objects.csv`
(output of `04_FindObjects`).

**Execution environment:** RSP / USDF (LSST kernel, Firefly)  
**README:** [`03_VisualiseCoadds/usdf_butler/README.md`](03_VisualiseCoadds/usdf_butler/README.md)

---

### `04_FindObjects/usdf_butler/`

**Topic:** Cross-match SIMBAD stable stars against the Rubin/LSSTCam
**object catalogue** (coadd-level photometry) via the Gen3 Butler.

**What it does:**
- Probes the Butler schema on a representative patch to discover which
  photometric columns are available, then loads object tables
  **patch-by-patch** (not tract-by-tract) to stay within the RSP memory
  limit.
- Cross-matches with `astropy SkyCoord.match_to_catalog_sky` (default
  radius 1.2 arcsec) and enriches each target with PSF fluxes, aperture
  fluxes, shape/morphology parameters, blendedness, and flags for all
  six bands.
- Converts PSF fluxes (nJy) to AB magnitudes; produces diagnostic plots:
  separation histogram, positional-offset scatter, match-status pie chart,
  r-band magnitude histogram, and a g−r vs r−i colour–colour diagram
  annotated with SIMBAD identifiers.
- Provides a preliminary `00_QueryButler.ipynb` scratchpad for exploring
  available dataset types, collections, and tract footprints before
  running the main pipeline.

**Key output:** `data_MATCHOBJ_02_out/targets_matched_lsst_objects.{csv,parquet}`
— consumed by `03_VisualiseCoadds` and `05_FindSources`.

**Execution environment:** RSP / USDF (LSST kernel)  
**README:** [`04_FindObjects/usdf_butler/README.md`](04_FindObjects/usdf_butler/README.md)

---

### `05_FindSources/usdf_butler/`

**Topic:** Extract per-visit **source-catalogue light curves** (`psfFlux`
vs. `expMidptMJD`) for stable stars with known spectral types, using the
Gen3 Butler.

**What it does:**
- Cross-matches each SIMBAD target against the LSSTCam `source` catalogue
  (single-visit detections), loaded memory-safely one patch at a time via
  the shared library `libExtractLightcurves.py` (incremental CSV writes,
  `gc.collect()` after each patch).
- Merges the resulting light curves with the DP2 visit table to attach
  `expMidptMJD`, airmass, seeing, and magnitude limit for each observation.
- Produces per-star **GridSpec 6 × 2 diagnostic figures** for each of the
  six bands: `psfFlux` vs. MJD time series (with Gaussian-fit mean and ±σ
  band overlaid) paired with a flux histogram and Gaussian fit.
- Includes a utility notebook (`99_ReadDP2VisitsTable.ipynb`) to cache the
  full DP2 visit table (with image-quality columns) as an ECSV file.

**Key outputs:** per-star light-curve Parquet files in
`data_MergeVisits_02_out/per_star/` and diagnostic figures in
`figs_PlotLC_03/`.

**Execution environment:** RSP / USDF (LSST kernel)  
**README:** [`05_FindSources/usdf_butler/README.md`](05_FindSources/usdf_butler/README.md)

---

## Data flow between series

```
SIMBAD API
    │
    ▼
01_askSimbad_stars
    │  data_SIMBAD_04/stable_stars_DDF<field>.csv
    │  data_SIMBAD_02/summary_visit_counts_per_star.parquet
    ▼
04_FindObjects/usdf_butler
    │  data_MATCHOBJ_02_out/targets_matched_lsst_objects.{csv,parquet}
    ├──────────────────────────────────────────────┐
    ▼                                              ▼
05_FindSources/usdf_butler             03_VisualiseCoadds/usdf_butler
    │  per-star psfFlux vs MJD             deepCoadd + star overlays
    │
    (uses also)
02_GetRubinVisits
    │  *_WithTracts.parquet visit tables
    │  (tract/patch + MJD + airmass + seeing)
```

---

## Common Butler configuration

| Parameter | `04_FindObjects` / `03_VisualiseCoadds` | `05_FindSources` | `02_GetRubinVisits` |
|-----------|----------------------------------------|-----------------|---------------------|
| Repository | `dp2_prep` | `dp2_prep` | `main` |
| Skymap | `lsst_cells_v2` | `lsst_cells_v2` | `lsst_cells_v2` |
| Collection | `LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881` | `LSSTCam/runs/DRP/v30_0_4/DM-54249` | `LSSTCam/defaults` |
| Instrument | LSSTCam | LSSTCam | LSSTCam |

---

## DDF fields

| Key | Field name | RA (deg) | Dec (deg) |
|-----|-----------|---------|----------|
| COSMOS | COSMOS Deep Drilling Field | 150.119 | +2.206 |
| XMM-LSS | XMM-LSS Deep Drilling Field | 35.708 | −4.750 |
| ELAIS-S1 | ELAIS-S1 Deep Drilling Field | 9.450 | −44.000 |
| ECDFS | Extended Chandra Deep Field South | 53.125 | −28.100 |
| EDFS-a | Euclid Deep Field South (pointing a) | 58.900 | −49.315 |
| EDFS-b | Euclid Deep Field South (pointing b) | 63.600 | −47.600 |
