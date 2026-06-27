# 04_FindObjects / usdf_butler

**Butler-based pipeline to identify and characterise photometrically stable
Simbad stars in Rubin/LSSTCam Deep Drilling Fields (DDFs)**

These notebooks run on the **Rubin Science Platform (RSP) at USDF** using the
Gen3 Butler and the DP2 data release
(`LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881`, skymap `lsst_cells_v2`).
They are designed to be executed **in numerical order**: each notebook writes
output files consumed by the next one.

---

## Notebooks

### `00_QueryButler.ipynb` — Butler registry exploration tool

An **interactive scratchpad** for inspecting the contents of a Butler
repository before running the main pipeline.

Key operations:
- Connect to the `dp2_prep` repository and list all registered dataset types.
- Filter dataset types by keyword (e.g. `object`, `table`, `deepCoadd`) to
  identify the correct names for downstream notebooks.
- Query available collections, instruments, and skymap configurations.
- Visualise tract/patch footprints on the sky using `lsst.sphgeom` and
  matplotlib polygon overlays.

> **Use this notebook first** whenever the collection or pipeline version
> changes, to verify which dataset-type names are available.

**Reference:** [104_2_Explore_the_butler_repo.ipynb](https://github.com/lsst/tutorial-notebooks/blob/main/DP1/100_How_to_Use_RSP_Tools/104_Butler_data_access/104_2_Explore_the_butler_repo.ipynb)

---

### `01_LoadLSSTCamDeepcoaddsCutout.ipynb` — DeepCoadd cutout viewer  (`NB_TAG = DEEPCCUTOUTS_01`)

Loads `deepCoadd` postage-stamp cutouts centred on Simbad targets and
displays them interactively in **Firefly** (the RSP image viewer).

Key operations:
- Reads the target list from `data_DEEPCCUTOUTS_01_in/summary_visit_counts_per_star_V17-21_r2.0deg.csv`.
- For each target, resolves the tract/patch with `skymap.findTract` /
  `findPatch` and retrieves the corresponding `deepCoadd` (or equivalent)
  dataset via the Butler.
- Applies `astropy` WCS and stretch normalisation (asinh, Z-scale) for
  matplotlib display.
- Sends cutouts to the Firefly RSP viewer via `lsst.afw.display` /
  `firefly_client`.
- Saves figures (PDF + PNG) to `figs_DEEPCCUTOUTS_01/`.

**Input directory:** `data_DEEPCCUTOUTS_01_in/`  
**Output directory:** `figs_DEEPCCUTOUTS_01/`

---

### `02_MatchTargetsWithLSSTCamObjects.ipynb` — Simbad × LSSTCam object cross-match  (`NB_TAG = MATCHOBJ_02`)

Cross-matches each Simbad target against the LSSTCam **object catalogue**
and enriches the output table with all available Rubin photometric
measurements (fluxes, flags, morphology).

Key operations:

1. **Schema inspection** — loads one representative patch to discover
   which `OBJ_COLUMNS` are actually present in the collection; missing
   columns are reported and skipped gracefully.
2. **Patch-level Butler loading** — object tables are loaded
   patch-by-patch (not tract-by-tract) and cached in memory, preventing
   RSP out-of-memory crashes (~1 M objects per tract vs ~5–20 k per patch).
3. **Nearest-neighbour cross-match** — `astropy.coordinates.SkyCoord.match_to_catalog_sky`
   with a configurable search radius (default 1.2 arcsec).
4. **Column enrichment** — every matched target receives all available
   `rubin_*` columns from `OBJ_COLUMNS`:
   - Global: `tract`, `patch`, `coord_ra/dec`, `coord_flag`, `ebv`
   - Per band (u, g, r, u, z, y): aperture fluxes (`ap12/17/25/50Flux`),
     PSF flux (`psfFlux`), per-band coordinates, `sizeExtendedness`,
     `blendedness`, `invalidPsfFlag`, `inputCount`
   - Unmatched rows receive `NaN` for all Rubin columns (uniform schema).
5. **AB magnitude conversion** — PSF fluxes (nJy) are converted to AB
   magnitudes using the nJy zero-point (31.4 mag) for bands g, r, i.
6. **Diagnostic plots** (saved as PDF + PNG):
   - Cross-match separation histogram
   - RA/Dec positional offset scatter plot
   - Match-status pie chart
   - PSF r-band magnitude histogram
   - **Colour-colour diagram g−r vs r−i**, colour-coded by r magnitude,
     with Simbad identifier annotations

**Input:** `data_DEEPCCUTOUTS_01_in/summary_visit_counts_per_star_V17-21_r2.0deg.csv`  
**Output directory:** `data_MATCHOBJ_02_out/`

| Output file | Content |
|-------------|---------|
| `targets_matched_lsst_objects.csv` / `.parquet` | All targets, all match statuses, full set of `rubin_*` columns |
| `targets_matched_only.csv` / `.parquet` | Matched-only subset for downstream notebooks |

**Reference:** [201_1_Object_table.ipynb](https://github.com/lsst/tutorial-notebooks/blob/main/DP1/200_Data_Products/201_Catalogs/201_1_Object_table.ipynb)

---

## Directory layout

```
usdf_butler/
├── README.md                          # this file
├── 00_QueryButler.ipynb               # Butler registry exploration
├── 01_LoadLSSTCamDeepcoaddsCutout.ipynb   # DeepCoadd cutout viewer
├── 02_MatchTargetsWithLSSTCamObjects.ipynb # Simbad × object cross-match
│
├── data_DEEPCCUTOUTS_01_in/           # Input CSV (shared by NB 01 and NB 02)
│   └── summary_visit_counts_per_star_V17-21_r2.0deg.csv
├── data_MATCHOBJ_02_out/              # Cross-match output tables
│   ├── targets_matched_lsst_objects.{csv,parquet}
│   └── targets_matched_only.{csv,parquet}
└── figs_DEEPCCUTOUTS_01/              # Figures from NB 01
```

---

## Butler configuration

| Parameter | Value |
|-----------|-------|
| Repository | `dp2_prep` |
| Collections | `LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage{1,2,3,4}` |
| Instrument | `LSSTCam` |
| Skymap | `lsst_cells_v2` |
| Pipeline | DP2 / Science Pipelines v30.0.0 |

---

## Dependencies

All notebooks run in the **LSST kernel** on the RSP (Python 3.12,
Science Pipelines stack). Key packages:

- `lsst.daf.butler`, `lsst.geom`, `lsst.sphgeom`, `lsst.afw.display`
- `lsst.skymap`
- `firefly_client` (NB 01 only)
- `astropy`, `numpy`, `pandas`, `matplotlib`

---

## Authors

- **Sylvie Dagoret-Campagne** — IJCLab/IN2P3/CNRS, Université Paris-Saclay
- Created: 2026-06-25 | Last update: 2026-06-27
