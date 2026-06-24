# 01_askSimbad_stars — Photometrically Stable Stars in Rubin/LSST Deep Drilling Fields

**Author:** Sylvie Dagoret-Campagne — IJCLab / IN2P3 / CNRS — Université Paris-Saclay  
**Created:** 2026-06-18  
**Last updated:** 2026-06-24  
**Context:** Rubin/LSST DP2 — Photometric calibration & systematic studies in Deep Drilling Fields

---

## Overview

This directory implements a **pipeline for identifying photometrically stable stars**
in the four primary Rubin/LSST Deep Drilling Fields (DDFs) for use as calibration
anchors in Rubin Butler Gen3 analyses.

The pipeline queries the [Simbad](https://simbad.cds.unistra.fr) astronomical
database for non-variable stars with known spectral types, cross-matches them
against the Rubin visit table to determine how many times each star has been
observed per band, retrieves optical sky image cutouts via the
[CDS HiPS2FITS](https://alasky.cds.unistra.fr/hips-image-services/hips2fits)
service for visual inspection, and finally exports per-field target lists in CSV
format for upload to the Rubin Science Platform (RSP) Portal.

### Deep Drilling Fields covered

| Field | RA (deg) | Dec (deg) |
|-------|----------|-----------|
| COSMOS | 150.1191 | +2.2058 |
| ELAIS-S1 | 9.4500 | −44.000 |
| XMM-LSS | 35.7080 | −4.750 |
| ECDFS | 53.1250 | −27.800 |

Search radius: **1.5 deg** from each field centre (adjustable).

### Magnitude selection

V-band magnitude range: **V ≤ 20** (bright enough for high-SNR calibration;
above Rubin saturation limit ~ 16 mag). The upper bound was updated from V ≤ 22
to V ≤ 20 on 2026-06-23 to focus on stars detectable with good SNR per visit.

---

## Notebook pipeline

The four main notebooks form a linear pipeline; each one consumes the output of
the previous one.

```
01_findStarsinSimbad.ipynb
        │
        └─► data_SIMBAD_01/master_stable_stars_*.csv
                │
                ▼
        02_StableStars_inRubinVisits.ipynb
                │
                └─► data_SIMBAD_02/summary_visit_counts_per_star.parquet
                        │
                        ├─► 03_show_StableStars_cutouts.ipynb  (visual QC)
                        │
                        └─► 04_makelist_forRubinPortal.ipynb
                                │
                                └─► data_SIMBAD_04/stable_stars_DDF<field>.csv
```

---

### `01_findStarsinSimbad.ipynb` — Simbad catalogue query

Queries the Simbad database for photometrically stable, non-variable stellar
objects in each DDF and builds a merged master catalogue.

**Simbad fields retrieved:**

| Field | Description |
|-------|-------------|
| `V`, `B`, `R`, `I` | Johnson magnitudes |
| `sp_type` | MK spectral type |
| `mesvar.vartyp` | Variability type from Simbad measurements table |
| `plx_value` | Parallax (mas) — Hipparcos / Gaia |
| `pmra`, `pmdec` | Proper motion components (mas/yr) |

**Selection criteria:**

- Object otype starts with `*` (stellar objects only — excludes galaxies, QSOs, etc.).
- Simbad variability flag is absent or blank (non-variable; see note below).
- V-band magnitude in the configured range (currently V ≤ 20).

> **Implementation note — masked-array pitfall:**  
> `mesvar.vartyp` is returned by Simbad as an astropy `MaskedColumn`.
> For stars with no variability record the entry is a **masked scalar**, not the
> literal string `'--'`. The correct test is `numpy.ma.is_masked()` first,
> not a string comparison against `'--'`, which would incorrectly flag all
> unobserved stars as variable.

**Spectral-type parsing:**  
Raw Simbad spectral types (e.g. `K2III`, `G8.5IV-V`, `M5V`, `DA2`) are parsed
into three components: MK temperature class (O–M), numeric subtype, and
luminosity class (I–V). Distributions are plotted per DDF.

**Diagnostic plots:**  
Sky distribution per DDF (RA/Dec scatter coloured by V mag), V-magnitude
histograms, MK class bar charts, luminosity class breakdown, and parallax
distributions (foreground / MW contamination check).

**Output files** (`data_SIMBAD_01/`):

| File | Description |
|------|-------------|
| `<field>_stars_V<min>-<max>_r<radius>deg.csv` | Per-field Simbad query cache |
| `master_stable_stars_V17-20_r1.5deg.csv` | Merged master catalogue (all DDFs) |

---

### `02_StableStars_inRubinVisits.ipynb` — Cross-match with Rubin visit table

For each star with a known spectral type from the master catalogue, determines
which Rubin visits covered it using a cone search on the visit pointing.

**Input visit table:**  
`rubindata_visits/visitTable-2025041500138-2026061900757_N88665_WithTracts.parquet`  
(~88 665 science visits, bands u/g/r/i/z/y, with `tract`, `patch`, `zenith_angle`).

**Cross-match strategy:**  
A star is matched to a visit if the angular separation between the star's
sky coordinates and the visit pointing (RA, Dec) is less than the Rubin
field-of-view radius (`MATCH_RADIUS_DEG = 1.75 deg`). The full visit
`SkyCoord` array is built once before the star loop for efficiency.

**Airmass computation:**  
Plane-parallel approximation: `X = 1 / cos(zenith_angle_rad)`, computed from
the `zenith_angle` column of the visit table. Airmass vs. MJD curves are
produced for each star.

**Outputs (`data_SIMBAD_02/`):**

| File / Directory | Description |
|------------------|-------------|
| `per_star/<simbad_id>.csv` | Per-star matched visit list with airmass |
| `per_star/<simbad_id>.parquet` | Same, in Parquet format |
| `summary_visit_counts_per_star.csv` | Summary table: one row per star, visit counts per band |
| `summary_visit_counts_per_star.parquet` | Same, in Parquet format |

The summary table columns include `simbad_id`, `spectral_type`, `mk_class_simple`,
`V_mag`, `ra_deg`, `dec_deg`, `field`, `n_visits_total`, and `n_u` … `n_y`
(per-band visit counts). Rows are sorted by `n_visits_total` descending.

**Figures** (`figs_SIMBAD_02/`): sky distribution of matched stars, per-star
visit-count bar charts, airmass vs. MJD time-series per star and per band.

---

### `03_show_StableStars_cutouts.ipynb` — HiPS2FITS sky image cutouts

For each stable star in the summary catalogue, retrieves a small sky image
cutout from the **CDS HiPS2FITS** web service and displays it for visual
quality control.

**HiPS surveys queried (PanSTARRS DR1):**

| HiPS ID | Band |
|---------|------|
| `CDS/P/PanSTARRS/DR1/g` | g |
| `CDS/P/PanSTARRS/DR1/r` | r |
| `CDS/P/PanSTARRS/DR1/i` | i |
| `CDS/P/PanSTARRS/DR1/z` | z |
| `CDS/P/PanSTARRS/DR1/y` | y |

**Cutout parameters (adjustable):**  
Field of view: 60 arcsec; image size: 512 × 512 pixels.

**Service details:**  
`astroquery.hips2fits.hips2fits.query()` with parameters `ra`, `dec`, `fov`,
`width`, `height`, `projection = TAN`. Falls back to a direct HTTP request
via `requests` if `astroquery.hips2fits` is unavailable.

**Display:** Each star is shown as a multi-panel figure (one panel per HiPS
band) with ZScale normalisation (`astropy.visualization.ZScaleInterval`) and
`AsinhStretch`. The star position is marked with a crosshair.

**Figures** (`figs_SIMBAD_03/`): multi-band cutout grids, one figure per star.

---

### `04_makelist_forRubinPortal.ipynb` — Export target lists for the RSP Portal

Prepares per-field CSV target lists suitable for upload to the
[Rubin Science Platform Portal](https://dp1.lsst.io/tutorials/portal/101/portal-101-2.html)
for cone-search queries and light-curve extraction.

**Selection applied:**  
Stars with a known spectral type **and** a minimum number of total Rubin visits
(`MIN_VISITS`, configurable). The summary table from notebook `02` is used as
the primary input; the NB01 master CSV is used as a fallback.

**Output format** (`data_SIMBAD_04/`):  
One CSV file per DDF with three columns: `simbad_id`, `ra`, `dec`.
Files are named `stable_stars_DDF<field>.csv` (e.g. `stable_stars_DDFCOSMOS.csv`).

These files can be loaded directly in the RSP Portal's
"Upload Table" / cone-search interface to retrieve Rubin
`diaObject` or `Object` photometry for each stable star.

---

### `99_testHips.ipynb` — HiPS2FITS development sandbox

Standalone notebook used to prototype and validate calls to the
`astroquery.hips2fits` interface before integrating them into notebook `03`.

Tests both the `query_with_wcs()` (full WCS specification) and the simpler
`query()` (RA, Dec, FoV, pixel size) interfaces against the COSMOS field centre.
Validates the list of available HiPS survey IDs (DSS2, PanSTARRS DR1, DES DR2).

This notebook does **not** depend on any previously produced data files.

---

## Data-flow summary

| Step | Input | Key operation | Output |
|------|-------|---------------|--------|
| **NB01** | Simbad API | Cone search per DDF, stability & magnitude filter | `data_SIMBAD_01/master_stable_stars_*.csv` |
| **NB02** | NB01 master CSV + Rubin visit Parquet | Angular cross-match, airmass computation | `data_SIMBAD_02/summary_visit_counts_per_star.parquet` |
| **NB03** | NB02 summary Parquet | CDS HiPS2FITS cutout retrieval | `figs_SIMBAD_03/` cutout grids |
| **NB04** | NB02 summary Parquet | Filter by visit count, export per-DDF CSVs | `data_SIMBAD_04/stable_stars_DDF<field>.csv` |

---

## Directory structure

| Path | Contents |
|------|----------|
| `data_SIMBAD_01/` | Per-field Simbad query caches + master stable-star catalogue |
| `data_SIMBAD_02/` | Per-star visit lists (CSV + Parquet) + summary table |
| `data_SIMBAD_04/` | Per-field portal target lists (CSV) |
| `figs_SIMBAD_02/` | Sky distribution, airmass curves, visit-count charts |
| `figs_SIMBAD_03/` | Multi-band HiPS cutout grids per star |
| `figs_SIMBAD_04/` | Diagnostic figures from the portal list notebook |
| `rubindata_visits/` | Rubin visit table Parquet (input, not version-controlled) |

---

## Software environment

| Item | Value |
|------|-------|
| Python | 3.13 |
| Kernel | `conda_py313` |
| Key packages | `astropy`, `astroquery`, `numpy`, `pandas`, `matplotlib`, `scipy` |
| Simbad access | `astroquery.simbad.Simbad` |
| HiPS2FITS access | `astroquery.hips2fits` (fallback: `requests`) |
| CDS HiPS2FITS URL | `https://alasky.cds.unistra.fr/hips-image-services/hips2fits` |
