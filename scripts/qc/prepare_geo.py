#!/usr/bin/env python3

"""Prepare tissue-level FASTQ symlinks and sample table for GEO submission."""

import argparse
import os
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd


TRACKING_REQUIRED = {"tissue", "original_rfid", "rfid", "status"}
FASTQS_REQUIRED = {"tissue", "rfid", "fastq"}


def validate_columns(path: Path, table: pd.DataFrame, required: set[str]) -> None:
    missing = required - set(table.columns)
    if missing:
        joined = ", ".join(sorted(missing))
        raise ValueError(f"{path} is missing required columns: {joined}")


def build_description(row: pd.Series) -> str:
    status = str(row["status"]).strip()
    original_rfid = str(row["original_rfid"]).strip()
    rfid = str(row["rfid"]).strip()

    if status != "included":
        return status
    if original_rfid != rfid:
        return f"Original label before genotype mismatch fix: {original_rfid}"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare GEO FASTQ links and sample sheet for one tissue"
    )
    parser.add_argument("tissue", help="Tissue ID matching the `tissue` column")
    parser.add_argument(
        "--tracking",
        type=Path,
        default=Path("samples/sample_tracking.tsv"),
        help="Tracking TSV (default: samples/sample_tracking.tsv)",
    )
    parser.add_argument(
        "--fastqs",
        type=Path,
        default=Path("samples/sample_fastqs.tsv"),
        help="FASTQ TSV (default: samples/sample_fastqs.tsv)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        help="Override FASTQ symlink directory (default: samples/geo_fastq-{tissue})",
    )
    parser.add_argument(
        "--out-samples",
        type=Path,
        help="Override GEO sample TSV path (default: samples/geo_samples-{tissue}.tsv)",
    )
    args = parser.parse_args()

    tissue = args.tissue.strip()
    if not tissue:
        raise ValueError("Tissue must be non-empty")

    tracking = pd.read_csv(args.tracking, sep="\t", dtype=str).fillna("")
    fastqs = pd.read_csv(args.fastqs, sep="\t", dtype=str).fillna("")
    validate_columns(args.tracking, tracking, TRACKING_REQUIRED)
    validate_columns(args.fastqs, fastqs, FASTQS_REQUIRED)

    tracking_tissue = tracking.loc[tracking["tissue"].str.strip() == tissue].copy()
    fastqs_tissue = fastqs.loc[fastqs["tissue"].str.strip() == tissue].copy()

    if tracking_tissue.empty:
        raise ValueError(f"No rows found for tissue '{tissue}' in {args.tracking}")
    if fastqs_tissue.empty:
        raise ValueError(f"No rows found for tissue '{tissue}' in {args.fastqs}")

    tracking_tissue["rfid"] = tracking_tissue["rfid"].astype(str).str.strip()
    fastqs_tissue["rfid"] = fastqs_tissue["rfid"].astype(str).str.strip()
    fastqs_tissue["fastq"] = fastqs_tissue["fastq"].astype(str).str.strip()
    if (tracking_tissue["rfid"] == "").any():
        raise ValueError(f"Empty RFIDs are not allowed in {args.tracking} for tissue '{tissue}'")
    if (fastqs_tissue["rfid"] == "").any():
        raise ValueError(f"Empty RFIDs are not allowed in {args.fastqs} for tissue '{tissue}'")

    geo_rows = tracking_tissue.copy()
    dup_rfids = geo_rows["rfid"][geo_rows["rfid"].duplicated()].unique()
    if len(dup_rfids) > 0:
        joined = ", ".join(sorted(dup_rfids))
        raise ValueError(f"Duplicate RFIDs in {args.tracking} for {tissue}: {joined}")

    tracked_rfids = set(geo_rows["rfid"].tolist())
    fastqs_geo = fastqs_tissue.copy()
    missing_fastqs = sorted(tracked_rfids - set(fastqs_geo["rfid"]))
    if missing_fastqs:
        joined = ", ".join(missing_fastqs)
        raise ValueError(f"Samples in {args.tracking} missing FASTQ entries in {args.fastqs}: {joined}")
    untracked_fastqs = sorted(set(fastqs_geo["rfid"]) - tracked_rfids)
    if untracked_fastqs:
        joined = ", ".join(untracked_fastqs)
        raise ValueError(
            f"FASTQ entries in {args.fastqs} are missing from {args.tracking}: {joined}"
        )

    # 1) Symlink each FASTQ in a flat directory, enforcing unique basenames.
    fastq_rel_paths = [x for x in fastqs_geo["fastq"] if x]
    basenames = [Path(p).name for p in fastq_rel_paths]
    duplicates = sorted(name for name, n in Counter(basenames).items() if n > 1)
    if duplicates:
        dup_str = ", ".join(duplicates)
        raise ValueError(
            "Cannot build flat GEO FASTQ directory because these file names repeat "
            f"across subdirectories: {dup_str}"
        )

    out_dir = args.out_dir or Path(f"samples/geo_fastq-{tissue}")
    out_dir.mkdir(parents=True, exist_ok=True)

    for rel_path in fastq_rel_paths:
        src = Path("v4") / tissue / "fastq" / rel_path
        if not src.exists():
            raise ValueError(f"Expected FASTQ is missing: {src}")
        dst = out_dir / Path(rel_path).name

        if dst.exists() or dst.is_symlink():
            if not dst.is_symlink():
                raise ValueError(f"Refusing to overwrite non-symlink destination: {dst}")
            current = dst.readlink()
            if current.is_absolute():
                current_abs = current.resolve()
            else:
                current_abs = (dst.parent / current).resolve()
            src_abs = src.resolve()
            if current_abs == src_abs:
                continue
            raise ValueError(
                f"Symlink already exists with different target: {dst} -> {current} (expected {src})"
            )

        src_for_link = Path(os.path.relpath(src, start=dst.parent))
        dst.symlink_to(src_for_link)

    # 2) Emit GEO sample table with rfid, description, fastq1..fastqN.
    geo_rows["description"] = geo_rows.apply(build_description, axis=1)

    files_by_rfid: dict[str, list[str]] = defaultdict(list)
    for _, row in fastqs_geo.iterrows():
        rfid = str(row["rfid"]).strip()
        fastq = str(row["fastq"]).strip()
        if rfid and fastq:
            files_by_rfid[rfid].append(Path(fastq).name)

    max_files = 0
    if not geo_rows.empty:
        max_files = max((len(files_by_rfid.get(rfid, [])) for rfid in geo_rows["rfid"]), default=0)

    out_records: list[dict[str, str]] = []
    file_cols = [f"fastq{i}" for i in range(1, max_files + 1)]
    for _, row in geo_rows.iterrows():
        rfid = row["rfid"]
        files = files_by_rfid.get(rfid, [])
        record: dict[str, str] = {"rfid": rfid, "description": row["description"]}
        for i, col in enumerate(file_cols):
            record[col] = files[i] if i < len(files) else ""
        out_records.append(record)

    out_table = pd.DataFrame(out_records, columns=["rfid", "description", *file_cols])
    out_samples = args.out_samples or Path(f"samples/geo_samples-{tissue}.tsv")
    out_samples.parent.mkdir(parents=True, exist_ok=True)
    out_table.to_csv(out_samples, sep="\t", index=False)


if __name__ == "__main__":
    main()
