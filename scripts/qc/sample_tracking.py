#!/usr/bin/env python3

"""Build a sample tracking table for rat IDs across datasets."""

import argparse
from pathlib import Path

import pandas as pd


def load_tissues(path: Path) -> list[str]:
    """Load dataset names from a one-column text file."""
    with path.open() as f:
        return [line.strip() for line in f if line.strip()]


def load_fastq_rfid_counts(path: Path) -> dict[str, int]:
    """Count fastq_map rows per RFID.

    Expected fastq_map formats per row:
    - single-end: fastq, rfid
    - paired-end: fastq1, fastq2, rfid
    """
    counts: dict[str, int] = {}
    with path.open() as f:
        for line_no, line in enumerate(f, start=1):
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) == 2:
                _, rfid = fields
            elif len(fields) == 3:
                _, _, rfid = fields
            else:
                raise ValueError(
                    f"{path}:{line_no} has {len(fields)} columns; expected 2 or 3"
                )
            counts[rfid] = counts.get(rfid, 0) + 1
    return counts


def load_rat_ids(path: Path) -> set[str]:
    """Load v4 release rat IDs from rat_ids.txt."""
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
    where updated_rfid may be empty for removed samples.
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


def make_removed_row(
    tissue: str,
    original_rfid: str,
    status: str,
    in_fastq_map: bool = False,
    n_fastq_rows: int = 0,
) -> dict[str, object]:
    """Create a canonical removed-sample row with no retained RFID."""
    return {
        "tissue": tissue,
        "original_rfid": original_rfid,
        "rfid": "",
        "in_fastq_map": in_fastq_map,
        "n_fastq_rows": n_fastq_rows,
        "is_genotyped": False,
        "in_v4_release": False,
        "status": status,
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Assemble a sample tracking table across datasets using "
            "{version}/{tissue}/fastq_map.txt and rat_ids.txt"
        )
    )
    parser.add_argument(
        "--version",
        type=str,
        default="v4",
        help="Version directory containing per-tissue inputs (default: v4)",
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
            "mixup relabeling; empty updated_rfid marks removed samples"
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
    for tissue in tissues:
        tissue_dir = Path(args.version) / tissue
        fastq_map = tissue_dir / "fastq_map.txt"
        rat_ids_file = tissue_dir / "rat_ids.txt"

        fastq_counts = load_fastq_rfid_counts(fastq_map)
        fastq_rfids = set(fastq_counts)
        v4_rfids = load_rat_ids(rat_ids_file)
        removed_originals = sorted(
            original
            for (relabel_tissue, original), updated in relabels.items()
            if relabel_tissue == tissue and not updated
        )
        relabel_targets = {
            updated
            for (relabel_tissue, _), updated in relabels.items()
            if relabel_tissue == tissue and updated
        }
        # Skip base rows for samples explicitly removed with no genotype match,
        # unless that RFID is now used as the kept label of another relabeled sample.
        skip_base_rfids = {
            original for original in removed_originals if original not in relabel_targets
        }
        all_rfids = sorted((fastq_rfids | v4_rfids) - skip_base_rfids)

        for rfid in all_rfids:
            original_rfid = updated_to_original.get((tissue, rfid), rfid)
            in_fastq_map = rfid in fastq_rfids
            in_v4_release = rfid in v4_rfids
            if in_fastq_map and in_v4_release:
                status = "included"
            elif in_fastq_map:
                status = "__missing_exclusion_reason__"
            else:
                status = "__missing_fastq_for_v4__"

            rows.append(
                {
                    "tissue": tissue,
                    "original_rfid": original_rfid,
                    "rfid": rfid,
                    "in_fastq_map": in_fastq_map,
                    "n_fastq_rows": fastq_counts.get(rfid, 0),
                    "is_genotyped": rfid in genotyped_rfids,
                    "in_v4_release": in_v4_release,
                    "status": status,
                }
            )

        # Keep explicit rows for samples removed due to no genotype match.
        for original_rfid in removed_originals:
            rows.append(
                make_removed_row(
                    tissue,
                    original_rfid,
                    "removed_no_genotype_match",
                    in_fastq_map=original_rfid in fastq_rfids,
                    n_fastq_rows=fastq_counts.get(original_rfid, 0),
                )
            )

    # Apply additional removal reasons from removed.tsv with ambiguity handling.
    by_tissue_original: dict[tuple[str, str], list[int]] = {}
    by_tissue_rfid: dict[tuple[str, str], list[int]] = {}
    for i, row in enumerate(rows):
        by_tissue_original.setdefault((row["tissue"], row["original_rfid"]), []).append(i)
        by_tissue_rfid.setdefault((row["tissue"], row["rfid"]), []).append(i)

    for (tissue, removed_rfid), reason in removed_reasons.items():
        status = f"removed_{reason}"
        orig_matches = by_tissue_original.get((tissue, removed_rfid), [])
        rfid_matches = by_tissue_rfid.get((tissue, removed_rfid), [])
        target_idx: int | None = None

        # Prefer exact original-sample matches, since removals apply to sample instances.
        if len(orig_matches) == 1:
            target_idx = orig_matches[0]
        elif len(orig_matches) > 1:
            raise ValueError(
                f"Ambiguous removed.tsv match for tissue={tissue} rfid={removed_rfid}: "
                "multiple original_rfid rows"
            )
        elif (tissue, removed_rfid) in updated_to_original:
            # RFID is currently used as a relabeled target for another sample.
            # Represent the removed original sample explicitly with blank final RFID.
            new_idx = len(rows)
            rows.append(make_removed_row(tissue, removed_rfid, status))
            by_tissue_original.setdefault((tissue, removed_rfid), []).append(new_idx)
            by_tissue_rfid.setdefault((tissue, ""), []).append(new_idx)
            continue
        elif len(rfid_matches) == 1:
            target_idx = rfid_matches[0]
        elif len(rfid_matches) > 1:
            raise ValueError(
                f"Ambiguous removed.tsv match for tissue={tissue} rfid={removed_rfid}: "
                "multiple rfid rows"
            )
        else:
            new_idx = len(rows)
            rows.append(make_removed_row(tissue, removed_rfid, status))
            by_tissue_original.setdefault((tissue, removed_rfid), []).append(new_idx)
            by_tissue_rfid.setdefault((tissue, ""), []).append(new_idx)
            continue

        rows[target_idx]["status"] = status
        if rows[target_idx]["in_v4_release"]:
            raise ValueError(
                "removed.tsv entry resolves to a v4-included sample: "
                f"tissue={tissue} rfid={removed_rfid} matched_row_rfid={rows[target_idx]['rfid']}"
            )
    for row in rows:
        if row["status"] == "__missing_exclusion_reason__":
            raise ValueError(
                "Sample in fastq_map but absent from v4 without explicit removal reason: "
                f"tissue={row['tissue']} original_rfid={row['original_rfid']} rfid={row['rfid']}"
            )
        if row["status"] == "__missing_fastq_for_v4__":
            raise ValueError(
                "Sample in v4 rat_ids but absent from fastq_map: "
                f"tissue={row['tissue']} original_rfid={row['original_rfid']} rfid={row['rfid']}"
            )

    table = pd.DataFrame(rows)
    table = table.sort_values(["tissue", "original_rfid", "rfid"]).reset_index(drop=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.out, sep="\t", index=False)


if __name__ == "__main__":
    main()
