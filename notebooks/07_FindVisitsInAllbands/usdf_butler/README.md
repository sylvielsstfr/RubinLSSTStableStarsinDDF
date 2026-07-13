# FindVisitsInAllbands / usdf_butler

Notebooks that find and characterize photometrically stable stars from
LSSTCamSources, per Deep Drilling Field (DDF) and per LSST band, using the
Rubin Science Platform Butler at USDF.

**Author:** Sylvie Dagoret-Campagne
**Affiliation:** IJCLab/IN2P3/CNRS, Universite Paris-Saclay

## Goal

For each DDF and each LSST band, find bright point-like sources (stars)
with `17 <= mag <= 22`, cross-match repeated detections of the same
physical object across visits (same RA/Dec within `MATCH_RADIUS_ARCSEC`),
keep objects with at least `MIN_VISITS_PER_BAND` visits in that band, and
compute the relative flux scatter `sigma_F/F` (in mmag) from `psfFlux`.
Processing is done band by band, DDF by DDF to control memory usage.

Expectation: the relative photometric scatter should be largest in **u**
and **y**, the bands most affected by atmospheric extinction variability
(Rayleigh + aerosols in u, water vapour/aerosols in y) and lower system
throughput.

## Notebooks

### 01_FindLSSTCamSourcesInAllbands.ipynb
Butler-facing notebook. Queries LSSTCamSources at USDF, cross-matches
detections into objects, computes per-object statistics
(`sigmaF_over_F_phot`, `sigmaF_over_F_meas`, `mmag_meas`, `mmag_phot`,
...), and writes one `objectstats_band_<band>.parquet` file per band to
`data_FindLSSTCamSources_01/`.

### 01b_FindLSSTCamSourcesInAllbands_keeptimeandfield.ipynb
Variant of 01 that additionally keeps the visit MJD and DDF field for each
detection, enabling light-curve-style diagnostics. Writes both
`objectstats_band_<band>.parquet` (per-object stats) and
`lcdeviation_band_<band>.parquet` (per-detection flux deviation vs. MJD)
to `data_FindLSSTCamSources_01b/`. 

### 02_ReadLSSTCamSourcesInAllbands.ipynb
Offline companion to 01. Does **not** touch the Butler: reads back the
`objectstats_band_<band>.parquet` files from `data_FindLSSTCamSources_01/`
and reproduces summary plots (boxplots, 2D histograms, per-magnitude-bin
histograms, error-bar and violin plots of `mmag_meas` vs. magnitude and
band, each compared against the photon-noise-only expectation
`mmag_phot`). Figures are written to `figs_PlotLSSTCamSources_02/`.

### 02b_ReadLSSTCamSourcesInAllbands.ipynb
Offline companion to 01b. Reads back the parquet files from
`data_FindLSSTCamSources_01b/` (both per-object stats and per-detection
light-curve deviations) and reproduces the same summary plots as 02, plus
light-curve diagnostics vs. MJD (per-visit median flux deviation and
dispersion, standard-error error bars, etc.). Figures are written to
`figs_ReadLSSTCamSourcesInAllbands_02b/`.

## Directory layout

```
usdf_butler/
├── 01_FindLSSTCamSourcesInAllbands.ipynb
├── 01b_FindLSSTCamSourcesInAllbands_keeptimeandfield.ipynb
├── 02_ReadLSSTCamSourcesInAllbands.ipynb
├── 02b_ReadLSSTCamSourcesInAllbands.ipynb
├── data_FindLSSTCamSources_01/       # parquet output of 01, input of 02
├── data_FindLSSTCamSources_01b/      # parquet output of 01b, input of 02b
├── figs_PlotLSSTCamSources_02/       # figures produced by 02
└── figs_ReadLSSTCamSourcesInAllbands_02b/  # figures produced by 02b
```

## Workflow

1. Run `01_...ipynb` (or `01b_...ipynb` for the time/field-aware variant)
   at USDF to query the Butler and write the per-band parquet files.
2. Run `02_...ipynb` (or `02b_...ipynb`) — locally or at USDF — to reload
   the parquet files and produce the diagnostic figures, without touching
   the Butler again.

## Key output columns (`objectstats_band_<band>.parquet`)

`object_id, n_visits, ra, dec, flux_mean, flux_std, mag_mean/mag_median,
sigmaF_over_F_phot, sigmaF_over_F_meas, mmag_meas, mmag_phot, ddf, band`

- `mmag_meas` — measured relative flux scatter across repeated visits.
- `mmag_phot` — photon-noise-only expectation for the same object; the
  floor that `mmag_meas` is compared against in the notebook-02 plots.
