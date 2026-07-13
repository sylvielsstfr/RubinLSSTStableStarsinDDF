#!/usr/bin/env python
"""
03_extract_lightcurves.py
=========================
Memory-efficient version of notebook 03_MatchTargetsWithLSSTCamSources.

Strategy to stay well below the RSP 16 GB limit:
  - NO in-memory cache of source tables: each (visit, detector) table is
    loaded, the single matching row is extracted, and the table is immediately
    released.
  - Results are written to disk incrementally (one CSV row per match) so the
    accumulator list stays tiny.
  - Matplotlib is imported only for the final diagnostic plots, which are
    generated from the already-saved CSV (small DataFrame).

Author: Sylvie Dagoret-Campagne
Affiliation: IJCLab/IN2P3/CNRS, Universite Paris-Saclay
Created: 2026-06-27
Updated: 2026-07-13
"""

import argparse
import gc
import logging
import os
import re

import astropy.units as u
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
from lsst.daf.butler import Butler, Timespan

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ===========================================================================
# CONFIGURATION  (edit this block, or override via CLI arguments)
# ===========================================================================

REPO = "dp2_prep"

COLLECTIONS = [
    "LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage1",
    "LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage2",
    "LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage3",
    "LSSTCam/runs/DRP/DP2/v30_0_0/DM-53881/stage4",
]

SKYMAP_NAME = "lsst_cells_v2"
INSTRUMENT = "LSSTCam"

NB_TAG = "MATCHSRC_01"
DIR_DATA_IN = "data_DEEPCCUTOUTS_01_in"
DIR_DATA_OUT = f"data_{NB_TAG}_out"
DIR_FIGS = f"figs_{NB_TAG}"

# TARGET_FILE = "summary_visit_counts_per_star_V17-21_r2.0deg.csv"
TARGET_FILE = "master_stable_stars_V17-21_r2.0deg.csv"

MATCH_RADIUS_ARCSEC = 1.0  # maximum angular separation for a valid match

DATE_START = "2025-04-01T00:00:00"
DATE_STOP = "2026-07-31T00:00:00"

# Photometric columns to extract from each source table.
# Keep this list as short as you need — every extra column costs RAM and I/O.
SRC_COLUMNS = [
    "coord_ra",
    "coord_dec",
    "parentSourceId",
    "x",
    "y",
    "xErr",
    "yErr",
    "ra",
    "dec",
    "raErr",
    "decErr",
    "calibFlux",
    "calibFluxErr",
    "ap09Flux",
    "ap09FluxErr",
    "ap09Flux_flag",
    "ap12Flux",
    "ap12FluxErr",
    "ap12Flux_flag",
    "ap17Flux",
    "ap17FluxErr",
    "ap17Flux_flag",
    "ap25Flux",
    "ap25FluxErr",
    "ap25Flux_flag",
    "ap35Flux",
    "ap35FluxErr",
    "ap35Flux_flag",
    "sky",
    "skyErr",
    "psfFlux",
    "psfFluxErr",
    "extendedness",
    "sizeExtendedness",
    "apFlux_12_0_flag",
    "apFlux_12_0_instFlux",
    "apFlux_12_0_instFluxErr",
    "apFlux_17_0_flag",
    "apFlux_17_0_instFlux",
    "apFlux_17_0_instFluxErr",
    "apFlux_35_0_flag",
    "apFlux_35_0_instFlux",
    "apFlux_35_0_instFluxErr",
    "extendedness_flag",
    "sizeExtendedness_flag",
    "localBackground_instFlux",
    "localBackground_instFluxErr",
    "localBackground_flag",
    "sky_source",
    "visit",
    "detector",
    "band",
    "physical_filter",
    "sourceId",
]


# ===========================================================================
# HELPERS
# ===========================================================================


def safe_name(s: str) -> str:
    """Replace characters that are invalid in filenames."""
    return re.sub(r"[^\w\-]", "_", s)


def find_col(df: pd.DataFrame, candidates: list[str]) -> str:
    """Return the first column name from *candidates* present in *df*."""
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(f"None of {candidates} found. Available columns (first 20): " f"{list(df.columns[:20])}")


def dataset_type_exists(butler: Butler, name: str) -> bool:
    """Return True if *name* is a registered dataset type in *butler*."""
    try:
        butler.registry.getDatasetType(name)
        return True
    except KeyError:
        return False


def probe_schema(
    butler: Butler,
    timespan: Timespan,
    ra: float,
    dec: float,
    src_columns: list[str],
) -> tuple[str, str, str, list[str]]:
    """Load the first available source table for position (ra, dec).

    Returns a tuple (ra_col, dec_col, id_col, src_columns_avail) where the
    column names are resolved from the actual schema and *src_columns_avail*
    is the subset of *src_columns* that actually exists in the table.
    """
    probe_refs = list(
        butler.query_datasets(
            "source",
            where=(
                "visit.timespan OVERLAPS :timespan AND "
                "visit_detector_region.region OVERLAPS POINT(:ra, :dec)"
            ),
            bind={"timespan": timespan, "ra": ra, "dec": dec},
        )
    )
    if not probe_refs:
        raise RuntimeError(
            f"No source refs found for the probe position "
            f"(ra={ra:.4f}, dec={dec:+.4f}).  "
            "Check the timespan and target coordinates."
        )

    log.info("Probe: %d refs found for first target", len(probe_refs))

    # Load only the columns we want — drastically reduces memory for the probe.
    df = butler.get(probe_refs[0], parameters={"columns": src_columns})
    if not isinstance(df, pd.DataFrame):
        df = df.to_pandas()

    ra_col = find_col(df, ["coord_ra", "ra", "RA", "ra_deg"])
    dec_col = find_col(df, ["coord_dec", "dec", "DEC", "dec_deg"])
    id_col = find_col(df, ["sourceId", "objectId", "object_id", "id"])
    avail = [c for c in src_columns if c in df.columns]

    log.info(
        "Schema probe  RA=%s  Dec=%s  ID=%s  avail_cols=%d/%d",
        ra_col,
        dec_col,
        id_col,
        len(avail),
        len(src_columns),
    )

    del df
    gc.collect()

    return ra_col, dec_col, id_col, avail


# ===========================================================================
# MAIN EXTRACTION LOOP
# ===========================================================================


def extract_lightcurves(
    butler: Butler,
    df_targets: pd.DataFrame,
    timespan: Timespan,
    ra_col: str,
    dec_col: str,
    id_col: str,
    src_columns_avail: list[str],
    match_radius_arcsec: float,
    out_csv: str,
    out_summary_csv: str,
) -> None:
    """Iterate over every target and every (visit, detector) ref.

    Performs the sky cross-match and writes matched rows to *out_csv*
    incrementally (append mode after each star) to keep RAM usage constant.

    Memory strategy
    ---------------
    - Source tables are loaded one at a time and released after the match.
    - No cache is kept in memory.
    - Matched rows are flushed to CSV after each star so the in-memory
      accumulator never grows beyond one star's worth of data.
    """
    # Columns already stored as dedicated metadata fields — skip when copying
    # photometric columns to avoid duplication.
    skip_in_row = {
        ra_col,
        dec_col,
        id_col,
        "visit",
        "detector",
        "band",
        "day_obs",
        "physical_filter",
    }
    photo_cols = [c for c in src_columns_avail if c not in skip_in_row]

    # Write CSV header (empty DataFrame with the right columns).
    header_row: dict = {
        "simbad_id": "",
        "target_ra": 0.0,
        "target_dec": 0.0,
        "visit": 0,
        "detector": 0,
        "band": "",
        "day_obs": 0,
        "physical_filter": "",
        "sep_arcsec": 0.0,
        "src_ra": 0.0,
        "src_dec": 0.0,
        "sourceId": 0,
    }
    for c in photo_cols:
        header_row[c] = 0.0
    pd.DataFrame([header_row]).head(0).to_csv(out_csv, index=False)

    summary_rows: list[dict] = []

    for idx, target in df_targets.iterrows():
        simbad_id = target["simbad_id"]
        ra_t = float(target["ra_deg"])
        dec_t = float(target["dec_deg"])
        # field = target["field"]

        # target coordinates
        tgt_sky = SkyCoord(ra=ra_t * u.deg, dec=dec_t * u.deg)

        log.info("[%3d] %s  ra=%.5f  dec=%+.5f", idx, simbad_id, ra_t, dec_t)

        # 1. Query refs
        try:
            refs = list(
                butler.query_datasets(
                    "source",
                    where=(
                        "visit.timespan OVERLAPS :timespan AND "
                        "visit_detector_region.region OVERLAPS POINT(:ra, :dec)"
                    ),
                    bind={"timespan": timespan, "ra": ra_t, "dec": dec_t},
                )
            )
        except Exception as exc:
            log.error("  ERROR querying refs: %s", exc)
            summary_rows.append(
                {
                    "simbad_id": simbad_id,
                    "ra_deg": ra_t,
                    "dec_deg": dec_t,
                    "n_refs": 0,
                    "n_matched": 0,
                    "n_failed": 0,
                    "status": "query_error",
                }
            )
            continue

        log.info("  -> %d refs (visit x detector pairs)", len(refs))

        n_matched = 0
        n_failed = 0
        batch: list[dict] = []

        for count, ref in enumerate(refs):
            did = ref.dataId
            visit = did["visit"]
            band = did.get("band", "?")
            day_obs = did.get("day_obs", -1)
            detector = did.get("detector", -1)
            physical_filter = did.get("physical_filter", "?")

            if count % 50 == 0:
                log.info(
                    "    ref %d/%d  visit=%s  band=%s  day_obs=%s",
                    count,
                    len(refs),
                    visit,
                    band,
                    day_obs,
                )

            # 2. Load source table (NO cache)
            df_src = None
            try:
                df_src = butler.get(ref, parameters={"columns": src_columns_avail})
                if not isinstance(df_src, pd.DataFrame):
                    df_src = df_src.to_pandas()
            except Exception as exc:
                log.warning(
                    "    WARNING: could not load ref (visit=%s det=%s): %s",
                    visit,
                    detector,
                    exc,
                )
                n_failed += 1
                continue

            if df_src is None or len(df_src) == 0:
                n_failed += 1
                del df_src
                continue

            # 3. SkyCoord catalogue for this detector
            ra_arr = df_src[ra_col].values
            dec_arr = df_src[dec_col].values

            # Heuristic: Butler stores coords in radians when max < 2pi + eps.
            unit_sky = u.rad if float(ra_arr.max()) <= 2 * np.pi + 0.1 else u.deg
            cat_sky = SkyCoord(ra=ra_arr * unit_sky, dec=dec_arr * unit_sky)

            # 4. Nearest-neighbour match
            best_i, sep2d, _ = tgt_sky.match_to_catalog_sky(cat_sky)
            sep_arcsec = float(sep2d.to(u.arcsec).value)

            if sep_arcsec > match_radius_arcsec:
                del df_src, ra_arr, dec_arr, cat_sky
                gc.collect()
                continue

            # 5. Extract the single matched row
            n_matched += 1
            matched = df_src.iloc[best_i]

            m_ra = float(matched[ra_col])
            m_dec = float(matched[dec_col])
            if unit_sky == u.rad:
                m_ra = np.degrees(m_ra)
                m_dec = np.degrees(m_dec)

            row: dict = {
                "simbad_id": simbad_id,
                "target_ra": ra_t,
                "target_dec": dec_t,
                "visit": visit,
                "detector": detector,
                "band": band,
                "day_obs": day_obs,
                "physical_filter": physical_filter,
                "sep_arcsec": sep_arcsec,
                "src_ra": m_ra,
                "src_dec": m_dec,
                "sourceId": (
                    int(matched[id_col]) if id_col in matched.index and pd.notna(matched[id_col]) else np.nan
                ),
            }
            for col in photo_cols:
                row[col] = matched.get(col, np.nan)

            batch.append(row)

            del df_src, ra_arr, dec_arr, cat_sky, matched
            gc.collect()

        # Flush this star's matches to disk
        if batch:
            pd.DataFrame(batch).to_csv(out_csv, mode="a", header=False, index=False)
            log.info("  flushed %d rows to disk", len(batch))
        batch.clear()

        log.info(
            "  matched: %d / %d refs  (load failures: %d)",
            n_matched,
            len(refs),
            n_failed,
        )

        summary_rows.append(
            {
                "simbad_id": simbad_id,
                "ra_deg": ra_t,
                "dec_deg": dec_t,
                "n_refs": len(refs),
                "n_matched": n_matched,
                "n_failed": n_failed,
                "status": "ok" if n_matched > 0 else "no_match",
            }
        )

    pd.DataFrame(summary_rows).to_csv(out_summary_csv, index=False)
    log.info("Match summary saved -> %s", out_summary_csv)


# ===========================================================================
# CHUNKED / STREAMING MEMORY-SAFE VARIANTS
#
# Rationale for the 16 GB USDF JupyterHub kernels:
#  1. `extract_lightcurves` above already streams matched rows to disk, but a
#     single kernel run still walks through *every* target and *every*
#     (visit, detector) ref with ONE Butler instance. The LSST Butler keeps
#     internal registry/datastore caches that grow across many thousands of
#     `.get()` calls and are NOT released by `gc.collect()` or `del`. Over a
#     long target list this is the real source of the OOM crash, not the
#     per-row Python objects (which are already freed).
#
#  2. The fix here is twofold:
#       a) let the notebook process only a CHUNK of the target list per
#          kernel run (`process_target_chunk`), so you can restart the
#          kernel between chunks and start from a clean process;
#       b) *within* a chunk, periodically recreate the Butler object
#          (`reset_butler_every`) so its internal caches are dropped even
#          without restarting the kernel.
#  3. Post-processing (CSV -> Parquet, per-star split, plots) is rewritten
#     to stream through the data in pandas `chunksize` pages instead of
#     loading the full light-curve table in memory, which was the *other*
#     major memory spike in the original notebook (`df_lc = pd.read_csv(...)`
#     followed by `df_lc.to_parquet(...)` and `save_per_star(df_lc, ...)`).
# ===========================================================================


def process_target_chunk(
    repo: str,
    collections: list[str],
    df_targets_chunk: pd.DataFrame,
    timespan: Timespan,
    ra_col: str,
    dec_col: str,
    id_col: str,
    src_columns_avail: list[str],
    match_radius_arcsec: float,
    out_csv: str,
    out_summary_csv: str,
    reset_butler_every: int = 5,
    resume: bool = True,
) -> None:
    """Memory-bounded, resumable variant of `extract_lightcurves`.

    Differences from `extract_lightcurves`:
      - Takes `repo`/`collections` instead of a pre-built `Butler` and
        creates/recreates the Butler itself, so its internal caches can be
        periodically released via `reset_butler_every`.
      - `out_csv`/`out_summary_csv` are expected to be CHUNK-SPECIFIC paths
        (e.g. suffixed with `_chunk00_of_10.csv`) so several kernel runs
        (one per chunk, ideally each in a fresh kernel) never collide.
      - If `resume=True` and `out_summary_csv` already exists (e.g. the
        previous run OOM'd partway through this chunk), targets already
        recorded there are skipped, and the summary is flushed to disk
        after *every* target (not just at the end) so a crash never loses
        more than one target's worth of progress.
    """
    skip_in_row = {
        ra_col,
        dec_col,
        id_col,
        "visit",
        "detector",
        "band",
        "day_obs",
        "physical_filter",
    }
    photo_cols = [c for c in src_columns_avail if c not in skip_in_row]

    # ── Resume bookkeeping ──────────────────────────────────────────────
    done_ids: set[str] = set()
    summary_rows: list[dict] = []
    if resume and os.path.exists(out_summary_csv):
        try:
            prev = pd.read_csv(out_summary_csv)
            summary_rows = prev.to_dict("records")
            done_ids = set(prev["simbad_id"].astype(str))
            log.info("Resume: %d targets already done in this chunk, skipping them", len(done_ids))
        except Exception as exc:
            log.warning("Could not read existing summary for resume (%s) - starting fresh", exc)

    # ── CSV header (only if the file does not already exist) ───────────
    if not os.path.exists(out_csv):
        header_row: dict = {
            "simbad_id": "",
            "target_ra": 0.0,
            "target_dec": 0.0,
            "visit": 0,
            "detector": 0,
            "band": "",
            "day_obs": 0,
            "physical_filter": "",
            "sep_arcsec": 0.0,
            "src_ra": 0.0,
            "src_dec": 0.0,
            "sourceId": 0,
        }
        for c in photo_cols:
            header_row[c] = 0.0
        pd.DataFrame([header_row]).head(0).to_csv(out_csv, index=False)

    butler = Butler(repo, collections=collections)
    log.info("Chunk processing: %d targets, reset_butler_every=%d", len(df_targets_chunk), reset_butler_every)

    n_processed_since_reset = 0
    for idx, target in df_targets_chunk.iterrows():
        simbad_id = target["simbad_id"]

        if resume and str(simbad_id) in done_ids:
            log.info("[%3d] %s  -> already processed, skipping (resume)", idx, simbad_id)
            continue

        if n_processed_since_reset > 0 and n_processed_since_reset % reset_butler_every == 0:
            log.info("  -- recreating Butler to release internal caches --")
            del butler
            gc.collect()
            butler = Butler(repo, collections=collections)

        ra_t = float(target["ra_deg"])
        dec_t = float(target["dec_deg"])
        tgt_sky = SkyCoord(ra=ra_t * u.deg, dec=dec_t * u.deg)

        log.info("[%3d] %s  ra=%.5f  dec=%+.5f", idx, simbad_id, ra_t, dec_t)

        try:
            refs = list(
                butler.query_datasets(
                    "source",
                    where=(
                        "visit.timespan OVERLAPS :timespan AND "
                        "visit_detector_region.region OVERLAPS POINT(:ra, :dec)"
                    ),
                    bind={"timespan": timespan, "ra": ra_t, "dec": dec_t},
                )
            )
        except Exception as exc:
            log.error("  ERROR querying refs: %s", exc)
            summary_rows.append(
                {
                    "simbad_id": simbad_id,
                    "ra_deg": ra_t,
                    "dec_deg": dec_t,
                    "n_refs": 0,
                    "n_matched": 0,
                    "n_failed": 0,
                    "status": "query_error",
                }
            )
            pd.DataFrame(summary_rows).to_csv(out_summary_csv, index=False)
            n_processed_since_reset += 1
            continue

        log.info("  -> %d refs (visit x detector pairs)", len(refs))

        n_matched = 0
        n_failed = 0
        batch: list[dict] = []

        for count, ref in enumerate(refs):
            did = ref.dataId
            visit = did["visit"]
            band = did.get("band", "?")
            day_obs = did.get("day_obs", -1)
            detector = did.get("detector", -1)
            physical_filter = did.get("physical_filter", "?")

            if count % 50 == 0:
                log.info(
                    "    ref %d/%d  visit=%s  band=%s  day_obs=%s", count, len(refs), visit, band, day_obs
                )

            df_src = None
            try:
                df_src = butler.get(ref, parameters={"columns": src_columns_avail})
                if not isinstance(df_src, pd.DataFrame):
                    df_src = df_src.to_pandas()
            except Exception as exc:
                log.warning("    WARNING: could not load ref (visit=%s det=%s): %s", visit, detector, exc)
                n_failed += 1
                continue

            if df_src is None or len(df_src) == 0:
                n_failed += 1
                del df_src
                continue

            ra_arr = df_src[ra_col].values
            dec_arr = df_src[dec_col].values
            unit_sky = u.rad if float(ra_arr.max()) <= 2 * np.pi + 0.1 else u.deg
            cat_sky = SkyCoord(ra=ra_arr * unit_sky, dec=dec_arr * unit_sky)

            best_i, sep2d, _ = tgt_sky.match_to_catalog_sky(cat_sky)
            sep_arcsec = float(sep2d.to(u.arcsec).value)

            if sep_arcsec > match_radius_arcsec:
                del df_src, ra_arr, dec_arr, cat_sky
                continue

            n_matched += 1
            matched = df_src.iloc[best_i]

            m_ra = float(matched[ra_col])
            m_dec = float(matched[dec_col])
            if unit_sky == u.rad:
                m_ra = np.degrees(m_ra)
                m_dec = np.degrees(m_dec)

            row: dict = {
                "simbad_id": simbad_id,
                "target_ra": ra_t,
                "target_dec": dec_t,
                "visit": visit,
                "detector": detector,
                "band": band,
                "day_obs": day_obs,
                "physical_filter": physical_filter,
                "sep_arcsec": sep_arcsec,
                "src_ra": m_ra,
                "src_dec": m_dec,
                "sourceId": (
                    int(matched[id_col]) if id_col in matched.index and pd.notna(matched[id_col]) else np.nan
                ),
            }
            for col in photo_cols:
                row[col] = matched.get(col, np.nan)

            batch.append(row)
            del df_src, ra_arr, dec_arr, cat_sky, matched

        if batch:
            pd.DataFrame(batch).to_csv(out_csv, mode="a", header=False, index=False)
            log.info("  flushed %d rows to disk", len(batch))
        batch.clear()
        gc.collect()

        log.info("  matched: %d / %d refs  (load failures: %d)", n_matched, len(refs), n_failed)

        summary_rows.append(
            {
                "simbad_id": simbad_id,
                "ra_deg": ra_t,
                "dec_deg": dec_t,
                "n_refs": len(refs),
                "n_matched": n_matched,
                "n_failed": n_failed,
                "status": "ok" if n_matched > 0 else "no_match",
            }
        )
        # Flush the summary after EVERY target so a crash never loses more
        # than the target currently in flight, and resume can pick up cleanly.
        pd.DataFrame(summary_rows).to_csv(out_summary_csv, index=False)
        n_processed_since_reset += 1

    del butler
    gc.collect()
    log.info("Chunk done. Match summary saved -> %s", out_summary_csv)


def merge_chunk_files(chunk_csv_paths: list[str], out_csv: str) -> None:
    """Concatenate several chunk light-curve CSVs into *out_csv* by streaming
    lines (no pandas involved), so RAM use stays independent of file size."""
    if not chunk_csv_paths:
        log.warning("merge_chunk_files: no chunk files given, nothing to do")
        return
    with open(out_csv, "w") as fout:
        for i, fpath in enumerate(chunk_csv_paths):
            with open(fpath) as fin:
                header = fin.readline()
                if i == 0:
                    fout.write(header)
                for line in fin:
                    fout.write(line)
    log.info("Merged %d chunk files -> %s", len(chunk_csv_paths), out_csv)


def merge_chunk_summaries(chunk_summary_paths: list[str], out_summary_csv: str) -> pd.DataFrame:
    """Concatenate the (small) per-chunk summary CSVs. These are tiny
    (one row per target) so loading them fully in memory is fine."""
    dfs = [pd.read_csv(p) for p in chunk_summary_paths if os.path.exists(p)]
    df_summary = pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()
    df_summary.to_csv(out_summary_csv, index=False)
    log.info("Merged %d chunk summaries -> %s (%d targets total)", len(dfs), out_summary_csv, len(df_summary))
    return df_summary


def csv_to_parquet_streaming(in_csv: str, out_parquet: str, chunksize: int = 200_000) -> None:
    """Convert a (potentially large) CSV to Parquet by reading it in pages,
    instead of `pd.read_csv(...).to_parquet(...)` which needs the whole
    table in RAM at once."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    writer = None
    n_rows = 0
    try:
        for chunk in pd.read_csv(in_csv, chunksize=chunksize):
            table = pa.Table.from_pandas(chunk, preserve_index=False)
            if writer is None:
                writer = pq.ParquetWriter(out_parquet, table.schema)
            writer.write_table(table)
            n_rows += len(chunk)
            del chunk, table
            gc.collect()
    finally:
        if writer is not None:
            writer.close()
    log.info("Streaming CSV -> Parquet done (%d rows) -> %s", n_rows, out_parquet)


def split_per_star_streaming(in_csv: str, dir_per_star: str, chunksize: int = 200_000) -> int:
    """Write one CSV per star, reading *in_csv* in pages instead of loading
    the full light-curve table in memory before grouping by star."""
    os.makedirs(dir_per_star, exist_ok=True)
    seen_stars: set[str] = set()
    for chunk in pd.read_csv(in_csv, chunksize=chunksize):
        for star_id, grp in chunk.groupby("simbad_id"):
            fname = safe_name(str(star_id))
            out_path = os.path.join(dir_per_star, f"{fname}_lc.csv")
            write_header = not os.path.exists(out_path)
            grp.to_csv(out_path, mode="a", header=write_header, index=False)
            seen_stars.add(star_id)
        del chunk
        gc.collect()
    log.info("Streaming per-star split done -> %d stars written to %s/", len(seen_stars), dir_per_star)
    return len(seen_stars)


def convert_per_star_to_parquet(dir_per_star: str) -> int:
    """Convert each per-star CSV (small - one star's light curve) to Parquet,
    one file at a time."""
    csv_files = sorted(f for f in os.listdir(dir_per_star) if f.endswith("_lc.csv"))
    for fname in csv_files:
        base = fname[:-4]
        df = pd.read_csv(os.path.join(dir_per_star, fname))
        df.to_parquet(os.path.join(dir_per_star, base + ".parquet"), index=False)
        del df
        gc.collect()
    log.info("Converted %d per-star CSV files to Parquet", len(csv_files))
    return len(csv_files)


def make_plots_streaming(
    in_csv: str,
    dir_per_star: str,
    dir_figs: str,
    match_radius: float,
    chunksize: int = 200_000,
    n_preview_stars: int = 6,
) -> None:
    """Diagnostic plots computed by streaming through *in_csv* in pages
    (histogram / band counts are accumulated incrementally) instead of
    requiring the full light-curve table in memory. The per-star flux
    preview reads only a handful of already-small per-star CSV files."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    os.makedirs(dir_figs, exist_ok=True)

    def savefig(name: str) -> None:
        for ext in ("pdf", "png"):
            plt.savefig(os.path.join(dir_figs, f"{name}.{ext}"), bbox_inches="tight")
        log.info("  -> saved %s.{pdf,png}", name)

    # ── Pass 1: accumulate separation histogram + per-star/band counts ──
    bin_edges = np.linspace(0, max(match_radius * 3, 1e-3), 31)
    hist_counts = np.zeros(len(bin_edges) - 1)
    band_counts: dict[tuple, int] = {}

    for chunk in pd.read_csv(in_csv, usecols=["simbad_id", "band", "sep_arcsec"], chunksize=chunksize):
        h, _ = np.histogram(chunk["sep_arcsec"].clip(upper=bin_edges[-1]), bins=bin_edges)
        hist_counts += h
        for key, n in chunk.groupby(["simbad_id", "band"]).size().items():
            band_counts[key] = band_counts.get(key, 0) + n
        del chunk
        gc.collect()

    fig, ax = plt.subplots(figsize=(6, 4))
    centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    ax.bar(centers, hist_counts, width=np.diff(bin_edges), edgecolor="k", linewidth=0.5)
    ax.axvline(match_radius, color="red", ls="--", label=f'search radius = {match_radius}"')
    ax.set_xlabel("Separation (arcsec)")
    ax.set_ylabel("Number of matches")
    ax.set_title("Cross-match separation - Simbad targets vs LSST sources (all visits)")
    ax.legend()
    plt.tight_layout()
    savefig("lc_crossmatch_separation_histogram")
    plt.close(fig)

    if band_counts:
        band_counts_series = pd.Series(band_counts)
        band_counts_df = band_counts_series.unstack(fill_value=0)
        fig, ax = plt.subplots(figsize=(max(6, len(band_counts_df) * 0.5), 4))
        band_counts_df.plot(kind="bar", ax=ax, width=0.8)
        ax.set_xlabel("Simbad target")
        ax.set_ylabel("Number of matched visits")
        ax.set_title("Visit count per star and band")
        ax.legend(title="band", bbox_to_anchor=(1.01, 1), loc="upper left")
        plt.xticks(rotation=45, ha="right", fontsize=7)
        plt.tight_layout()
        savefig("lc_visits_per_star_band")
        plt.close(fig)

    # ── psfFlux preview: only read the first N already-small per-star files ──
    bands_to_plot = ["r", "g", "i", "u", "z", "y"]
    band_colors = {"u": "purple", "g": "blue", "r": "green", "i": "orange", "z": "red", "y": "brown"}

    star_files = sorted(f for f in os.listdir(dir_per_star) if f.endswith("_lc.csv"))[:n_preview_stars]
    if not star_files:
        return

    dfs = [pd.read_csv(os.path.join(dir_per_star, f)) for f in star_files]
    flux_col = None
    for df in dfs:
        if "psfFlux" in df.columns:
            flux_col = "psfFlux"
            break
    if flux_col is None:
        for df in dfs:
            if "calibFlux" in df.columns:
                flux_col = "calibFlux"
                break

    if flux_col:
        n_stars = len(dfs)
        fig, axes = plt.subplots(n_stars, 1, figsize=(10, 3 * n_stars), sharex=False)
        if n_stars == 1:
            axes = [axes]

        for ax, df_star in zip(axes, dfs, strict=False):
            star_id = df_star["simbad_id"].iloc[0] if "simbad_id" in df_star.columns else "?"
            df_star = df_star.sort_values(["band", "visit"])
            for band in bands_to_plot:
                df_b = df_star[df_star["band"] == band]
                if len(df_b) == 0:
                    continue
                x = np.arange(len(df_b))
                y = df_b[flux_col].values
                yerr = (
                    df_b[flux_col + "Err"].values if flux_col + "Err" in df_b.columns else np.zeros(len(df_b))
                )
                ax.errorbar(
                    x,
                    y,
                    yerr=yerr,
                    fmt="o",
                    ms=3,
                    lw=0.8,
                    color=band_colors.get(band, "gray"),
                    label=f"{band} ({len(df_b)} pts)",
                )
            ax.set_title(star_id, fontsize=8)
            ax.set_xlabel("Visit index (proxy - MJD to be added)")
            ax.set_ylabel(f"{flux_col} [nJy]")
            ax.legend(fontsize=7, ncol=3)

        plt.suptitle(f"Light curves - {flux_col} (all bands)", y=1.01)
        plt.tight_layout()
        savefig("lc_psfFlux_preview")
        plt.close(fig)


# ===========================================================================
# POST-PROCESSING: save per-star files + diagnostic plots
# ===========================================================================


def save_per_star(df_lc: pd.DataFrame, dir_per_star: str) -> None:
    """Write one CSV + Parquet file per star into *dir_per_star*."""
    os.makedirs(dir_per_star, exist_ok=True)
    n = 0
    for star_id, grp in df_lc.groupby("simbad_id"):
        fname = safe_name(star_id)
        base = os.path.join(dir_per_star, f"{fname}_lc")
        grp.to_csv(base + ".csv", index=False)
        grp.to_parquet(base + ".parquet", index=False)
        n += 1
    log.info("Saved %d per-star LC files -> %s/", n, dir_per_star)


def make_plots(
    df_lc: pd.DataFrame,
    df_summary: pd.DataFrame,
    dir_figs: str,
    match_radius: float,
) -> None:
    """Generate and save diagnostic plots to *dir_figs*."""
    import matplotlib

    matplotlib.use("Agg")  # non-interactive backend (safe on RSP)
    import matplotlib.pyplot as plt

    os.makedirs(dir_figs, exist_ok=True)

    def savefig(name: str) -> None:
        """Save current figure as PDF and PNG."""
        for ext in ("pdf", "png"):
            plt.savefig(os.path.join(dir_figs, f"{name}.{ext}"), bbox_inches="tight")
        log.info("  -> saved %s.{pdf,png}", name)

    # Separation histogram
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(df_lc["sep_arcsec"], bins=30, edgecolor="k", linewidth=0.5)
    ax.axvline(match_radius, color="red", ls="--", label=f'search radius = {match_radius}"')
    ax.set_xlabel("Separation (arcsec)")
    ax.set_ylabel("Number of matches")
    ax.set_title("Cross-match separation - Simbad targets vs LSST sources (all visits)")
    ax.legend()
    plt.tight_layout()
    savefig("lc_crossmatch_separation_histogram")
    plt.close(fig)

    # Visit count per star/band
    band_counts = df_lc.groupby(["simbad_id", "band"]).size().unstack(fill_value=0)
    fig, ax = plt.subplots(figsize=(max(6, len(band_counts) * 0.5), 4))
    band_counts.plot(kind="bar", ax=ax, width=0.8)
    ax.set_xlabel("Simbad target")
    ax.set_ylabel("Number of matched visits")
    ax.set_title("Visit count per star and band")
    ax.legend(title="band", bbox_to_anchor=(1.01, 1), loc="upper left")
    plt.xticks(rotation=45, ha="right", fontsize=7)
    plt.tight_layout()
    savefig("lc_visits_per_star_band")
    plt.close(fig)

    # psfFlux preview (first 6 stars, all bands)
    bands_to_plot = ["r", "g", "i", "u", "z", "y"]
    band_colors = {
        "u": "purple",
        "g": "blue",
        "r": "green",
        "i": "orange",
        "z": "red",
        "y": "brown",
    }
    flux_col = (
        "psfFlux" if "psfFlux" in df_lc.columns else "calibFlux" if "calibFlux" in df_lc.columns else None
    )

    if flux_col:
        n_stars = min(6, df_lc["simbad_id"].nunique())
        stars = df_lc["simbad_id"].unique()[:n_stars]
        fig, axes = plt.subplots(n_stars, 1, figsize=(10, 3 * n_stars), sharex=False)
        if n_stars == 1:
            axes = [axes]

        for ax, star_id in zip(axes, stars, strict=False):
            df_star = df_lc[df_lc["simbad_id"] == star_id].sort_values(["band", "visit"])
            for band in bands_to_plot:
                df_b = df_star[df_star["band"] == band]
                if len(df_b) == 0:
                    continue
                x = np.arange(len(df_b))
                y = df_b[flux_col].values
                yerr = (
                    df_b[flux_col + "Err"].values if flux_col + "Err" in df_b.columns else np.zeros(len(df_b))
                )
                ax.errorbar(
                    x,
                    y,
                    yerr=yerr,
                    fmt="o",
                    ms=3,
                    lw=0.8,
                    color=band_colors.get(band, "gray"),
                    label=f"{band} ({len(df_b)} pts)",
                )
            ax.set_title(star_id, fontsize=8)
            ax.set_xlabel("Visit index (proxy - MJD to be added)")
            ax.set_ylabel(f"{flux_col} [nJy]")
            ax.legend(fontsize=7, ncol=3)

        plt.suptitle(f"Light curves - {flux_col} (all bands)", y=1.01)
        plt.tight_layout()
        savefig("lc_psfFlux_preview")
        plt.close(fig)


# ===========================================================================
# ENTRY POINT
# ===========================================================================


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    p = argparse.ArgumentParser(
        description=(
            "Extract per-visit source light curves from the Rubin Butler "
            "(memory-efficient, streaming output)."
        )
    )
    p.add_argument("--repo", default=REPO, help=f"Butler repository alias or path (default: {REPO})")
    p.add_argument(
        "--target-file",
        default=TARGET_FILE,
        help=f"Input CSV file name inside DIR_DATA_IN (default: {TARGET_FILE})",
    )
    p.add_argument(
        "--dir-data-in",
        default=DIR_DATA_IN,
        help=f"Input data directory (default: {DIR_DATA_IN})",
    )
    p.add_argument(
        "--dir-data-out",
        default=DIR_DATA_OUT,
        help=f"Output data directory (default: {DIR_DATA_OUT})",
    )
    p.add_argument(
        "--dir-figs",
        default=DIR_FIGS,
        help=f"Output figures directory (default: {DIR_FIGS})",
    )
    p.add_argument(
        "--match-radius",
        type=float,
        default=MATCH_RADIUS_ARCSEC,
        help=f"Cross-match radius in arcsec (default: {MATCH_RADIUS_ARCSEC})",
    )
    p.add_argument(
        "--date-start",
        default=DATE_START,
        help=f"Start date ISO8601 UTC (default: {DATE_START})",
    )
    p.add_argument(
        "--date-stop",
        default=DATE_STOP,
        help=f"Stop date ISO8601 UTC (default: {DATE_STOP})",
    )
    p.add_argument("--no-plots", action="store_true", help="Skip diagnostic plot generation")
    return p.parse_args()


def main() -> None:
    """Run the full light-curve extraction pipeline."""
    args = parse_args()

    os.makedirs(args.dir_data_out, exist_ok=True)
    os.makedirs(args.dir_figs, exist_ok=True)

    # Load targets
    target_path = os.path.join(args.dir_data_in, args.target_file)
    df_targets = pd.read_csv(target_path)
    log.info("Loaded %d targets from %s", len(df_targets), target_path)

    # Butler
    butler = Butler(args.repo, collections=COLLECTIONS)
    log.info("Butler initialised | repo: %s", args.repo)

    # Timespan (Butler requires TAI)
    t1 = Time(args.date_start, format="isot", scale="utc")
    t2 = Time(args.date_stop, format="isot", scale="utc")
    timespan = Timespan(
        Time(t1.mjd, format="mjd", scale="tai"),
        Time(t2.mjd, format="mjd", scale="tai"),
    )
    log.info("Timespan MJD [%.1f, %.1f]  (delta=%.0f days)", t1.mjd, t2.mjd, t2.mjd - t1.mjd)

    # Schema probe
    first = df_targets.iloc[0]
    ra_col, dec_col, id_col, src_cols_avail = probe_schema(
        butler,
        timespan,
        ra=float(first["ra_deg"]),
        dec=float(first["dec_deg"]),
        src_columns=SRC_COLUMNS,
    )

    # Output paths
    out_lc_csv = os.path.join(args.dir_data_out, "all_stars_lightcurves.csv")
    out_sum_csv = os.path.join(args.dir_data_out, "lightcurve_match_summary.csv")
    dir_per_star = os.path.join(args.dir_data_out, "per_star")

    # Main loop
    extract_lightcurves(
        butler=butler,
        df_targets=df_targets,
        timespan=timespan,
        ra_col=ra_col,
        dec_col=dec_col,
        id_col=id_col,
        src_columns_avail=src_cols_avail,
        match_radius_arcsec=args.match_radius,
        out_csv=out_lc_csv,
        out_summary_csv=out_sum_csv,
    )

    # Read back the full LC table (on disk -> small RAM footprint)
    df_lc = pd.read_csv(out_lc_csv)
    df_summary = pd.read_csv(out_sum_csv)
    log.info("df_lc: %d rows x %d cols", len(df_lc), len(df_lc.columns))

    # Parquet copy of the global LC table
    out_lc_parquet = os.path.join(args.dir_data_out, "all_stars_lightcurves.parquet")
    df_lc.to_parquet(out_lc_parquet, index=False)
    log.info("Parquet copy saved -> %s", out_lc_parquet)

    # Per-star files
    save_per_star(df_lc, dir_per_star)

    # Diagnostic plots
    if not args.no_plots and len(df_lc) > 0:
        make_plots(df_lc, df_summary, args.dir_figs, args.match_radius)

    log.info("Done.")


if __name__ == "__main__":
    main()
