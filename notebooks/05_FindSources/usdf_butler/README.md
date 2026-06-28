# 05_FindSources — Butler-based Light Curve Extraction for Stable Stars

**Author:** Sylvie Dagoret-Campagne  
**Affiliation:** IJCLab/IN2P3/CNRS, Université Paris-Saclay  
**Project:** Rubin/LSST — Photometrically Stable Stars in Deep Drilling Fields (DDFs)

---

## Overview

This notebook series extracts and analyses per-visit `source`-catalogue light
curves for a curated set of **photometrically stable stars selected from
SIMBAD for which the spectral type is known**.  Knowing the spectral type
allows the flux time series to be interpreted in a physically motivated way
(e.g. predicting the expected DCR shift, computing colour corrections, or
anchoring the atmospheric calibration).

> **Scope of this series.**  The notebooks documented here deal exclusively
> with SIMBAD targets whose spectral type is available.  A companion series
> (to be developed) will handle SIMBAD stable stars for which **no spectral
> type** is known; that series may rely on photometric typing or colour-based
> selection instead.

All computations run at the **USDF** (US Data Facility) using the Rubin Gen3
Butler and the DP2 data release
(`LSSTCam/runs/DRP/v30_0_4/DM-54249`, skymap `lsst_cells_v2`).

---

## Notebook Descriptions

### `01_MatchTargetsWithLSSTCamSources.ipynb`

**Goal:** Cross-match the SIMBAD stable-star target list against the LSSTCam
`source` catalogue retrieved via the Butler, and save per-star light-curve
files.

**Workflow:**

1. Reads the input target CSV (columns: `simbad_id`, `ra_deg`, `dec_deg`,
   spectral-type metadata) from `data_DEEPCCUTOUTS_01_in/`.
2. For each target, locates the Butler tract/patch covering the target
   coordinates.
3. Loads the corresponding `source` dataset (auto-discovered from the registry)
   using memory-safe incremental access (no full-table in-memory cache) via
   `libExtractLightcurves.py`.
4. Performs a nearest-neighbour sky cross-match with
   `astropy.coordinates.SkyCoord.match_to_catalog_sky` within a configurable
   search radius.
5. Retains a pruned set of columns (flux, flux error, shape, flags, `visitId`,
   `band`, …) to stay within the RSP 16 GB memory limit.
6. Writes results incrementally to CSV and writes final per-star Parquet files
   under `data_MATCHSRC_01_out/per_star/`.
7. Produces diagnostic figures in `figs_MATCHSRC_01/`.

**Key inputs:** `data_DEEPCCUTOUTS_01_in/` (target list)  
**Key outputs:** `data_MATCHSRC_01_out/` (global CSV + per-star CSV/Parquet)

---

### `02_MergeLCsourceswithMJD.ipynb`

**Goal:** Enrich the light-curve files produced by notebook `01` with precise
observation timestamps (`expMidptMJD`) and observing-condition metadata from
the official DP2 visit table.

**Workflow:**

1. Loads the DP2 visit table (`dp2_visits_table_with_iq.ecsv`, produced by
   `99_ReadDP2VisitsTable.ipynb`) and builds a lookup dictionary keyed on
   `visitId`.
2. For every light-curve file in `data_MATCHSRC_01_out/` (global and per-star),
   merges on the `visit` column to attach:
   - `expMidptMJD` — observation mid-point in Modified Julian Date,
   - `obsStartMJD` — shutter-open time,
   - `band_visit` — photometric band from the visit table (renamed to avoid
     collision with the source-catalogue `band` column),
   - `airmass`, `mean_seeing`, `mean_maglim`.
3. Writes the enriched files to `data_MergeVisits_02_out/per_star/`, mirroring
   the same sub-directory layout.
4. Saves summary figures in `figs_MergeVisits_02/`.

**Key inputs:** `data_MATCHSRC_01_out/`, `data_MergeVisits_02_in/` (visit table)  
**Key outputs:** `data_MergeVisits_02_out/` (time-stamped per-star light curves)

---

### `03_PlotLCwithMJD.ipynb`

**Goal:** Produce publication-quality diagnostic plots of `psfFlux` vs.
`expMidptMJD` for each stable star, one multi-panel figure per star.

**Workflow:**

For each per-star file in `data_MergeVisits_02_out/per_star/`, generates a
**GridSpec 6 × 2 figure**:

- **Left column (6 rows, one per band u g r i z y):**  
  `psfFlux` vs. `expMidptMJD` with individual error bars (`psfFluxErr`).  
  A shared x-axis range is enforced across all bands.  
  A secondary top x-axis shows the calendar date (YYYY-MM-DD).  
  The y-axis is clipped to ±N·σ_IQR around the median.  
  The Gaussian-fit mean is overlaid as a dashed horizontal line; the ±σ band
  is shown as a shaded region.

- **Right column (6 rows, one per band):**  
  Histogram of `psfFlux` for the same band, with a Gaussian fit overlaid;
  fit mean (μ) and standard deviation (σ) are annotated.

Bands with fewer than 3 good measurements are skipped (empty subplot).  
Figures are saved as both PDF and PNG in `figs_PlotLC_03/`.

**Key inputs:** `data_MergeVisits_02_out/per_star/`  
**Key outputs:** `figs_PlotLC_03/` (PDF + PNG diagnostic figures)

---

### `99_ReadDP2VisitsTable.ipynb`

**Goal:** Utility notebook — retrieve and cache the DP2 visit table for use
by the rest of the series.

Reads the consolidated DP2 visit table (including image-quality columns:
airmass, seeing, magnitude limit) from the Butler or a pre-cached ECSV file,
and writes it to `data_MergeVisits_02_in/dp2_visits_table_with_iq.ecsv`.
This notebook is run once as a prerequisite before running notebooks `01`–`03`.

**Key outputs:** `data_MergeVisits_02_in/dp2_visits_table_with_iq.ecsv`

---

## Support Library

### `libExtractLightcurves.py`

A Python module shared across the notebook series, providing memory-safe
helper functions for Butler-based source extraction:

| Function | Description |
|---|---|
| `safe_name(s)` | Sanitise a SIMBAD identifier into a filesystem-safe string. |
| `find_col(df, candidates)` | Return the first matching column name from a priority list. |
| `dataset_type_exists(butler, name, **kwargs)` | Probe the Butler registry without raising on missing datasets. |
| `probe_schema(butler, dataset_type, data_id)` | Inspect the column schema of a Butler dataset without loading the full table. |
| `save_per_star(df, out_dir, star_id)` | Incrementally write a per-star CSV row and flush to disk immediately (avoids accumulating large DataFrames in memory). |

The library is designed to stay well within the RSP 16 GB memory limit by
loading one `(visit, detector)` source table at a time, extracting the single
matching row, and immediately releasing the table with `gc.collect()`.

---

## Directory Layout

```
usdf_butler/
├── README.md                          # This file
├── libExtractLightcurves.py           # Shared memory-safe Butler helpers
│
├── 99_ReadDP2VisitsTable.ipynb        # Prerequisite: cache the DP2 visit table
├── 01_MatchTargetsWithLSSTCamSources.ipynb   # Cross-match targets → sources
├── 02_MergeLCsourceswithMJD.ipynb     # Attach MJD + observing conditions
├── 03_PlotLCwithMJD.ipynb             # Diagnostic light-curve figures
│
├── data_DEEPCCUTOUTS_01_in/           # Input: SIMBAD target list (with spectral types)
├── data_MATCHSRC_01_out/              # Output of notebook 01
├── data_MergeVisits_02_in/            # Input for notebook 02 (DP2 visit table)
├── data_MergeVisits_02_out/           # Output of notebook 02
│
├── figs_MATCHSRC_01/                  # Figures from notebook 01
├── figs_MergeVisits_02/               # Figures from notebook 02
└── figs_PlotLC_03/                    # Figures from notebook 03
```

---

## Recommended Execution Order

```
99_ReadDP2VisitsTable   →   01_MatchTargets   →   02_MergeLC   →   03_PlotLC
```

---

## Dependencies

- **Rubin Science Platform (RSP) / USDF** with the `lsst_distrib` stack
- Python ≥ 3.11 (tested with the `conda_py313` kernel)
- `astropy`, `numpy`, `pandas`, `matplotlib`, `scipy`
- Butler collection: `LSSTCam/runs/DRP/v30_0_4/DM-54249`
- Skymap: `lsst_cells_v2`

---

## Related Series

| Series | Description |
|---|---|
| `04_FindObjects/` | Butler-based object-catalogue cross-matching (coadd-level photometry). |
| `05_FindSources/` (this series) | Per-visit source-catalogue light curves — **SIMBAD stars with known spectral type**. |
| *(forthcoming)* | Per-visit light curves for SIMBAD stable stars **without known spectral type**. |
