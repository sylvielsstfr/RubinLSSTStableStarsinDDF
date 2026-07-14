# 01_MatchTargetsWithLSSTCamSources

Cross-matches Simbad targets against LSSTCam `source` catalogues (DP2, via the
USDF Butler) to build per-target light curves, then splits the result into
per-star files and produces diagnostic plots.

## Files

- `01_MatchTargetsWithLSSTCamSources.ipynb` — the notebook (config, schema
  probe, extraction loop, post-processing).
- `libExtractLightcurves.py` — all the logic used by the notebook:
  schema probing, the extraction loop, and post-processing helpers.
- `data_DEEPCCUTOUTS_01_in/` — input target list.
- `data_MATCHSRC_01_out/` — output CSV/Parquet light curves, per-star files,
  and figures.

## Memory-bounded design

Running on a 16 GB USDF JupyterHub kernel, the naive version of this notebook
runs out of memory. Two independent causes, and how each is addressed:

1. **Extraction loop.** The LSST `Butler` accumulates internal
   registry/datastore caches across many thousands of `.get()` calls; these
   are not released by `gc.collect()` or `del`. `process_target_chunk(...)`
   works around this by:
   - recreating the `Butler` every `reset_butler_every` targets (default 5),
   - flushing matched rows to disk after every target instead of
     accumulating them in memory,
   - flushing the summary CSV after every target and supporting
     `resume=True`, so if the kernel still crashes you just **restart the
     kernel and re-run the same cell** — already-completed targets (found in
     `out_sum_csv`) are skipped automatically.

2. **Post-processing.** The original approach (`pd.read_csv` of the full
   light-curve table, then `.to_parquet()` / `.groupby()` / plotting on it)
   loads everything in RAM at once. The notebook now uses streaming
   equivalents instead:
   - `csv_to_parquet_streaming` — CSV → Parquet via `pyarrow.ParquetWriter`,
     page by page (`chunksize=200_000`).
   - `split_per_star_streaming` — writes one CSV per star by reading the
     combined CSV in pages and appending each page's rows to the right
     per-star file.
   - `convert_per_star_to_parquet` — converts each (small) per-star CSV to
     Parquet one file at a time.
   - `make_plots_streaming` — separation histogram and per-star/band counts
     are accumulated incrementally across pages; the flux preview plot only
     reads a handful of the already-small per-star files.

The original non-streaming functions (`extract_lightcurves`, `save_per_star`,
`make_plots`) are still in `libExtractLightcurves.py` for reference / small
target lists, but the notebook no longer uses them.

## Usage

1. Run the config + schema-probe cells once.
2. Run the extraction cell (`process_target_chunk(...)`). If the kernel dies,
   just restart it and re-run this same cell — it resumes automatically.
3. Run the post-processing cells (streaming Parquet conversion, per-star
   split, diagnostic plots).

## Tuning

If the extraction loop still runs out of memory, lower `reset_butler_every`
(e.g. 2–3 instead of 5) to recreate the Butler more often, at the cost of
some extra re-initialization time.
