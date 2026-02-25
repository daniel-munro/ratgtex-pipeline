#!/usr/bin/env python3

"""Build a sample tracking table for rat IDs across datasets."""

import argparse
from pathlib import Path

import pandas as pd


def load_tissues(path: Path) -> list[str]:
    """Load dataset names from a one-column text file."""
    with path.open() as f:
        return [line.strip() for line in f if line.strip()]


def load_fastq_entries(path: Path) -> list[tuple[str, str]]:
    """Load one FASTQ file per row from fastq_map.

    Returns a list of (fastq, rfid). For paired-end rows, emits one row per file.
    """
    entries: list[tuple[str, str]] = []
    with path.open() as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) == 2:
                fastq, rfid = fields
                entries.append((fastq, rfid))
            elif len(fields) == 3:
                fastq1, fastq2, rfid = fields
                entries.append((fastq1, rfid))
                entries.append((fastq2, rfid))
            else:
                raise ValueError(
                    f"{path}:{line_no} has {len(fields)} columns; expected 2 or 3"
                )
    return entries


def load_rat_ids(path: Path) -> set[str]:
    """Load release rat IDs from a rat_ids file."""
    with path.open() as f:
        return {line.strip() for line in f if line.strip()}


def load_genotyped_rfids(path: Path) -> set[str]:
    """Load genotyped rat IDs from genotyping log CSV column `rfid`."""
    table = pd.read_csv(path, usecols=["rfid"], dtype=str)
    rfids = table["rfid"].dropna().astype(str).str.strip()
    return set(rfids.loc[rfids != ""])


def load_relabels(path: Path) -> dict[tuple[str, str], str]:
    """Load per-tissue sample relabeling decisions.

    Returns a mapping:
      (tissue, original_rfid) -> updated_rfid
    """
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"tissue", "original_rfid", "updated_rfid"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing))}"
        )

    mapping: dict[tuple[str, str], str] = {}
    for i, row in table.iterrows():
        tissue = row["tissue"].strip()
        original = row["original_rfid"].strip()
        updated = row["updated_rfid"].strip()
        if not tissue or not original:
            raise ValueError(
                f"{path}: row {i + 2} has empty tissue/original_rfid, which is invalid"
            )
        if not updated:
            raise ValueError(
                f"{path}: row {i + 2} has empty updated_rfid; "
                "use an explicit unknown-{tissue}-{n} label instead"
            )
        key = (tissue, original)
        if key in mapping and mapping[key] != updated:
            raise ValueError(
                f"{path}: conflicting relabel entries for tissue={tissue} original_rfid={original}"
            )
        mapping[key] = updated
    return mapping


def build_updated_to_original(
    relabels: dict[tuple[str, str], str],
) -> dict[tuple[str, str], str]:
    """Build reverse mapping for relabeled (non-empty updated_rfid) entries.

    Returns:
      (tissue, updated_rfid) -> original_rfid
    """
    reverse: dict[tuple[str, str], str] = {}
    for (tissue, original), updated in relabels.items():
        if not updated:
            continue
        key = (tissue, updated)
        if key in reverse and reverse[key] != original:
            raise ValueError(
                f"Conflicting originals map to tissue={tissue} updated_rfid={updated}"
            )
        reverse[key] = original
    return reverse


def load_removed(path: Path) -> dict[tuple[str, str], str]:
    """Load per-tissue removals with free-text reason.

    Returns:
      (tissue, rfid) -> reason
    """
    table = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    required = {"tissue", "rfid", "reason"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(
            f"{path} is missing required columns: {', '.join(sorted(missing))}"
        )

    mapping: dict[tuple[str, str], str] = {}
    for i, row in table.iterrows():
        tissue = row["tissue"].strip()
        rfid = row["rfid"].strip()
        reason = row["reason"].strip()
        if not tissue or not rfid:
            raise ValueError(
                f"{path}: row {i + 2} has empty tissue/rfid, which is invalid"
            )
        key = (tissue, rfid)
        if key in mapping and mapping[key] != reason:
            raise ValueError(
                f"{path}: conflicting removal entries for tissue={tissue} rfid={rfid}"
            )
        mapping[key] = reason
    return mapping


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Assemble a sample tracking table across datasets using "
            "samples/fastq_map_current/fastq_map-{tissue}.txt and "
            "samples/rat_ids_current/rat_ids-{tissue}.txt"
        )
    )
    parser.add_argument(
        "--tissues",
        type=Path,
        default=Path("tissues.dup.txt"),
        help="File listing datasets/tissues to include, one per line",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("samples/sample_tracking.tsv"),
        help="Output TSV file path",
    )
    parser.add_argument(
        "--out-fastqs",
        type=Path,
        default=Path("samples/sample_fastqs.tsv"),
        help="Output TSV file with columns tissue, rfid, fastq (one FASTQ per row)",
    )
    parser.add_argument(
        "--genotyping-log",
        type=Path,
        default=Path("geno/genotyping_log.csv"),
        help="CSV file containing genotyping metadata with column `rfid`",
    )
    parser.add_argument(
        "--relabeled",
        type=Path,
        default=Path("samples/relabeled.tsv"),
        help=(
            "TSV with columns tissue, original_rfid, updated_rfid for sample "
            "mixup relabeling; updated_rfid must be non-empty"
        ),
    )
    parser.add_argument(
        "--removed",
        type=Path,
        default=Path("samples/removed.tsv"),
        help=(
            "TSV with columns tissue, rfid, reason for additional removed samples; "
            "rfid is post-relabel sample ID"
        ),
    )
    args = parser.parse_args()

    tissues = load_tissues(args.tissues)
    genotyped_rfids = load_genotyped_rfids(args.genotyping_log)
    relabels = load_relabels(args.relabeled)
    updated_to_original = build_updated_to_original(relabels)
    removed_reasons = load_removed(args.removed)
    rows = []
    fastq_rows = []
    for tissue in tissues:
        fastq_map = Path("samples/fastq_map_current") / f"fastq_map-{tissue}.txt"
        rat_ids_file = Path("samples/rat_ids_current") / f"rat_ids-{tissue}.txt"

        fastq_entries = load_fastq_entries(fastq_map)
        fastq_counts: dict[str, int] = {}
        for fastq_file, raw_rfid in fastq_entries:
            fastq_counts[raw_rfid] = fastq_counts.get(raw_rfid, 0) + 1
            fastq_rows.append({"tissue": tissue, "rfid": raw_rfid, "fastq": fastq_file})
        fastq_rfids = set(fastq_counts)
        v4_rfids = load_rat_ids(rat_ids_file)
        missing_fastq_for_v4 = sorted(v4_rfids - fastq_rfids)
        if missing_fastq_for_v4:
            joined = ", ".join(missing_fastq_for_v4)
            raise ValueError(
                f"Samples in rat_ids-{tissue}.txt are missing from fastq_map-{tissue}.txt: {joined}"
            )
        for (relabel_tissue, original), updated in relabels.items():
            if relabel_tissue != tissue:
                continue
            if updated not in fastq_rfids:
                raise ValueError(
                    f"Relabel target missing from fastq_map-{tissue}.txt: "
                    f"original_rfid={original} updated_rfid={updated}"
                )

        all_rfids = sorted(fastq_rfids)
        for rfid in all_rfids:
            original_rfid = updated_to_original.get((tissue, rfid), rfid)
            in_v4_release = rfid in v4_rfids
            is_genotyped = rfid in genotyped_rfids
            reason = removed_reasons.get((tissue, rfid), "")
            if in_v4_release:
                if reason:
                    raise ValueError(
                        "Sample is in v4 but also listed in removed.tsv: "
                        f"tissue={tissue} rfid={rfid} reason={reason}"
                    )
                if not is_genotyped:
                    raise ValueError(
                        "Sample is in v4 but is not genotyped: "
                        f"tissue={tissue} rfid={rfid}"
                    )
                status = "included"
            elif reason:
                status = f"removed_{reason}"
            elif rfid.startswith("unknown-"):
                status = "removed_no_matching_genotype_sample"
            elif rfid.startswith("dup-"):
                status = "removed_duplicate"
            elif not is_genotyped:
                status = "removed_not_genotyped"
            else:
                raise ValueError(
                    "Sample absent from v4 is genotyped but has no explicit removal reason: "
                    f"tissue={tissue} original_rfid={original_rfid} rfid={rfid}"
                )
            rows.append(
                {
                    "tissue": tissue,
                    "original_rfid": original_rfid,
                    "rfid": rfid,
                    "n_fastq_rows": fastq_counts.get(rfid, 0),
                    "is_genotyped": is_genotyped,
                    "in_v4_release": in_v4_release,
                    "status": status,
                }
            )

    known_rfids = {(str(row["tissue"]), str(row["rfid"])) for row in rows}
    for tissue, removed_rfid in removed_reasons:
        if (tissue, removed_rfid) not in known_rfids:
            raise ValueError(
                "removed.tsv entry is not present in the fastq map for this tissue: "
                f"tissue={tissue} rfid={removed_rfid}"
            )

    table = pd.DataFrame(rows)
    table = table.sort_values(["tissue", "original_rfid", "rfid"]).reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)

    fastq_table = pd.DataFrame(fastq_rows, columns=["tissue", "rfid", "fastq"])
    fastq_table = fastq_table.sort_values(["tissue", "rfid", "fastq"]).reset_index(
        drop=True
    )
    args.out_fastqs.parent.mkdir(parents=True, exist_ok=True)
    fastq_table.to_csv(args.out_fastqs, sep="\t", index=False)


if __name__ == "__main__":
    main()
