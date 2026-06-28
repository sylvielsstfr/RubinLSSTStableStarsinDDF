# RubinLSSTStableStarsinDDF

**Identification and characterisation of photometrically stable stars
in Rubin/LSST Deep Drilling Fields for atmospheric and photometric calibration**

[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)
[![PyPI](https://img.shields.io/pypi/v/rubinlsststablestarsinddf?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/rubinlsststablestarsinddf/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/LSSTDESC/rubinlsststablestarsinddf/smoke-test.yml)](https://github.com/LSSTDESC/rubinlsststablestarsinddf/actions/workflows/smoke-test.yml)
[![Codecov](https://codecov.io/gh/LSSTDESC/rubinlsststablestarsinddf/branch/main/graph/badge.svg)](https://codecov.io/gh/LSSTDESC/rubinlsststablestarsinddf)
[![Read The Docs](https://img.shields.io/readthedocs/rubinlsststablestarsinddf)](https://rubinlsststablestarsinddf.readthedocs.io/)
[![Benchmarks](https://img.shields.io/github/actions/workflow/status/LSSTDESC/rubinlsststablestarsinddf/asv-main.yml?label=benchmarks)](https://LSSTDESC.github.io/rubinlsststablestarsinddf/)

**Authors:** Sylvie Dagoret-Campagne, Martin Rodriguez Monroy  
**Affiliation:** IJCLab / IN2P3 / CNRS — Université Paris-Saclay  
**GitHub:** [LSSTDESC/rubinlsststablestarsinddf](https://github.com/LSSTDESC/rubinlsststablestarsinddf)

---

## Scientific objectives

The Rubin Observatory Legacy Survey of Space and Time (LSST) will observe the
southern sky repeatedly over ten years. In the **Deep Drilling Fields** (DDFs),
hundreds of visits per band are accumulated, enabling high-precision time-domain
photometry. Reaching the sub-percent photometric calibration required for Type
Ia supernova cosmology and weak lensing demands a thorough understanding of two
dominant sources of systematic uncertainty:

1. **Atmospheric transmission variability** — precipitable water vapour (PWV),
   aerosol optical depth, and grey extinction fluctuate on timescales of minutes
   to seasons and vary across the focal plane.

2. **Differential Chromatic Refraction (DCR)** — the atmosphere refracts light
   differentially with wavelength, shifting the apparent position of every
   source along the parallactic angle direction by an amount that depends on
   colour and airmass. In difference imaging, a colour mismatch between the
   science image and the template produces dipole artefacts aligned with the
   zenith direction.

Photometrically **stable stars with known spectral types** are ideal calibration
anchors to measure and correct both effects, because:

- their intrinsic flux is constant and their spectral energy distribution is
  known from the MK classification,
- repeated Rubin observations give a long baseline light curve from which any
  residual variability can be attributed to the atmosphere or the instrument,
- their known colour predicts the expected DCR shift as a function of airmass
  and parallactic angle, providing a direct handle on the atmospheric refractivity.

This package implements a complete end-to-end pipeline that, starting from the
SIMBAD astronomical database, produces per-star multi-band photometric light
curves and deepCoadd visualisations for all six Rubin DDFs, ready for
atmospheric calibration analysis.

---

## Pipeline structure

The analysis is organised as five sequential notebook series located in
`notebooks/`. Each series produces output files consumed by the next.

```
SIMBAD (astroquery)
        │
        ▼
01_askSimbad_stars      Query SIMBAD for non-variable stars with known
        │               spectral types in each DDF. Build the master
        │               target catalogue with MK classification,
        │               Gaia parallaxes, V magnitudes, and per-band
        │               Rubin visit counts.
        │
        ▼
02_GetRubinVisits       Build LSSTCam visit tables enriched with
        │               Butler skymap tract/patch indices, MJD
        │               timestamps, airmass, seeing, and limiting
        │               magnitude (Butler registry and/or ConsDB).
        │
        ▼
04_FindObjects          Cross-match SIMBAD targets against the Rubin
        │               LSSTCam object catalogue (deepCoadd-level) via
        │               the Gen3 Butler. Attach PSF and aperture fluxes
        │               for all six bands (u g r i z y).
        │
        ├──────────────────────────────────────────┐
        ▼                                          ▼
05_FindSources                          03_VisualiseCoadds
Per-visit source-catalogue              Display deepCoadd images in
light curves: psfFlux vs.              Firefly with stable-star
expMidptMJD per band, with             overlays annotated with
Gaussian-fit diagnostics.              SIMBAD name and spectral type.
```

---

## Notebook series

| Directory | Topic | Environment |
|-----------|-------|-------------|
| [`notebooks/01_askSimbad_stars/`](notebooks/01_askSimbad_stars/README.md) | SIMBAD query, master star catalogue, visit cross-match | local Mac (`conda_py313`) |
| [`notebooks/02_GetRubinVisits/`](notebooks/02_GetRubinVisits/README.md) | LSSTCam visit table with tract/patch indices | RSP / USDF |
| [`notebooks/03_VisualiseCoadds/usdf_butler/`](notebooks/03_VisualiseCoadds/usdf_butler/README.md) | DeepCoadd visualisation with star overlays in Firefly | RSP / USDF |
| [`notebooks/04_FindObjects/usdf_butler/`](notebooks/04_FindObjects/usdf_butler/README.md) | Simbad × LSSTCam object catalogue cross-match | RSP / USDF |
| [`notebooks/05_FindSources/usdf_butler/`](notebooks/05_FindSources/usdf_butler/README.md) | Per-visit source light curve extraction | RSP / USDF |

A detailed description of each series, including the data flow between them,
is given in [`notebooks/README.md`](notebooks/README.md).

---

## Deep Drilling Fields

| Key | Field | RA (deg) | Dec (deg) |
|-----|-------|----------|-----------|
| COSMOS | COSMOS Deep Drilling Field | 150.119 | +2.206 |
| XMM-LSS | XMM-LSS Deep Drilling Field | 35.708 | −4.750 |
| ELAIS-S1 | ELAIS-S1 Deep Drilling Field | 9.450 | −44.000 |
| ECDFS | Extended Chandra Deep Field South | 53.125 | −28.100 |
| EDFS-a | Euclid Deep Field South (pointing a) | 58.900 | −49.315 |
| EDFS-b | Euclid Deep Field South (pointing b) | 63.600 | −47.600 |

---

## Butler configuration (RSP / USDF)

| Series | Repository | Collection | Skymap |
|--------|-----------|------------|--------|
| 03, 04 | `dp2_prep` | `LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage{1–4}` | `lsst_cells_v2` |
| 05 | `dp2_prep` | `LSSTCam/runs/DRP/v30_0_4/DM-54249` | `lsst_cells_v2` |
| 02 | `main` | `LSSTCam/defaults` | `lsst_cells_v2` |

---

## Repository layout

```
RubinLSSTStableStarsinDDF/
├── README.md                    # this file
├── pyproject.toml               # package metadata (LINCC-Frameworks template)
├── requirements.txt
├── LICENSE                      # MIT
│
├── src/rubinlsststablestarsinddf/   # Python package (in development)
│   ├── __init__.py
│   ├── example_module.py
│   └── example_benchmarks.py
│
├── notebooks/                   # Analysis notebook series (main deliverable)
│   ├── README.md                # Pipeline overview and data-flow diagram
│   ├── 01_askSimbad_stars/
│   ├── 02_GetRubinVisits/
│   ├── 03_VisualiseCoadds/
│   ├── 04_FindObjects/
│   └── 05_FindSources/
│
├── docs/                        # Sphinx / ReadTheDocs documentation
├── tests/                       # pytest test suite
└── benchmarks/                  # asv performance benchmarks
```

---

## Installation

### On the RSP / USDF (LSST kernel)

The LSST Science Pipelines stack is pre-installed. Clone the repository and
use it directly from the `notebooks/` directory — no additional installation
is required for the notebooks.

```bash
git clone https://github.com/LSSTDESC/rubinlsststablestarsinddf.git
```

### On a local machine (Mac / Linux)

```bash
conda create -n rubinstars python=3.11
conda activate rubinstars
git clone https://github.com/LSSTDESC/rubinlsststablestarsinddf.git
cd rubinlsststablestarsinddf
./.setup_dev.sh       # installs dev dependencies and pre-commit hooks
conda install pandoc  # optional, for Sphinx notebook rendering
```

The local notebooks (`01_askSimbad_stars`) run with the `conda_py313` kernel
and require: `astropy`, `astroquery`, `numpy`, `pandas`, `matplotlib`,
`scipy`, `healpy`.

---

## Related work


- **AuxTel / Spectractor PWV analysis** — seasonal and intra-night PWV
  measurements at Cerro Pachón using the LATISS spectrograph, providing the
  atmospheric water-vapour time series against which the stable-star flux
  residuals can be calibrated.


---

## License

MIT — see [`LICENSE`](LICENSE).
