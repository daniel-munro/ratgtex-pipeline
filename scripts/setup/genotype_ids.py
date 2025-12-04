#!/usr/bin/env python3
"""
Select genotyping sample IDs for a set of requested RFIDs.

Inputs:
  - A text file with one RFID per line (may contain underscores).
  - A text file with one VCF sample ID per line. Each VCF ID encodes an RFID in one of:
      ddGBS_{info1}_{info2}_{RFID}
      Riptide{info1}_{RFID}_{info2}
      Rattaca{info1}_{RFID}_{info2}
      LSWrnaseq_{RFID}_LSWNNNNN
    Where info1 and info2 never contain underscores, but RFID may contain underscores.

Outputs:
  - Filtered RFIDs: subset of requested RFIDs that have at least one VCF sample.
  - Selected VCF IDs: exactly one VCF ID per output RFID, using:
      - If there are Riptide matches for the RFID, choose the one with the highest integer
        in the first Riptide match group.
      - Else if there are ddGBS(R) matches, choose the last occurring ddGBS(R) match.
      - Else if there is an LSWrnaseq match, choose that one (at most one per RFID).
      - Else exclude the RFID (ignore Rattaca).

Usage:
  python scripts/setup/genotype_ids.py \
    --rfids geno/original/rfids_for_geno_11_1.txt \
    --vcf-ids round11_1_vcf_ids.txt \
    --out-tsv selected_rfids.tsv
"""

from __future__ import annotations

import argparse
import re
import sys
from typing import Dict, Iterable, List, Optional, Tuple


# Precompile regex patterns for speed and clarity.
_RE_DDGBS = re.compile(r"^ddGBSR?_([^_]+)_([^_]+)_(.+)$")
_RE_RIPTIDE = re.compile(r"^Riptide([^_]+)_(.+)_([^_]+)$")
_RE_LSWRNASEQ = re.compile(r"^LSWrnaseq_(.+)_LSWNNNNN$")


def select_vcf_per_rfid(
    requested_rfids: Iterable[str],
    vcf_ids: Iterable[str],
) -> Tuple[List[str], List[str]]:
    """Select one VCF ID per requested RFID.

    Selection policy:
      - Prefer Riptide matches, choosing the one with the highest integer from the
        first Riptide match group.
      - If no Riptide matches for an RFID, use the last occurring ddGBS(R) match.
      - If neither, exclude the RFID (ignore Rattaca).
    """
    req_list = list(requested_rfids)
    vcf_list = list(vcf_ids)

    # Track best Riptide (by highest rank) and last ddGBS per RFID.
    best_riptide_for_rfid: Dict[str, Tuple[int, str]] = {}  # rfid -> (rank, sid)
    last_ddgbs_for_rfid: Dict[str, str] = {}                # rfid -> sid
    lsw_rnaseq_for_rfid: Dict[str, str] = {}                # rfid -> sid
    unknown_patterns = 0

    for sid in vcf_list:
        s = sid.strip()

        # Try Riptide first, since it carries a rank.
        m_rip = _RE_RIPTIDE.match(s)
        if m_rip:
            # m_rip.group(1) is expected to be an integer per spec.
            try:
                rank = int(m_rip.group(1))
            except ValueError:
                unknown_patterns += 1
                continue
            rfid = m_rip.group(2)
            prev = best_riptide_for_rfid.get(rfid)
            if prev is None or rank > prev[0]:
                best_riptide_for_rfid[rfid] = (rank, sid)
            continue

        # Then ddGBS(R)
        m_dd = _RE_DDGBS.match(s)
        if m_dd:
            rfid = m_dd.group(3)
            last_ddgbs_for_rfid[rfid] = sid
            continue

        # LSW rnaseq fallback (only one expected per RFID).
        m_lsw = _RE_LSWRNASEQ.match(s)
        if m_lsw:
            rfid = m_lsw.group(1)
            lsw_rnaseq_for_rfid[rfid] = sid
            continue

        # Ignore Rattaca and anything else.
        unknown_patterns += 1

    filtered_rfids: List[str] = []
    selected_vcfs: List[str] = []
    riptide_selected = 0
    ddgbs_selected = 0
    lsw_selected = 0

    for rfid in req_list:
        if rfid in best_riptide_for_rfid:
            filtered_rfids.append(rfid)
            selected_vcfs.append(best_riptide_for_rfid[rfid][1])
            riptide_selected += 1
        elif rfid in last_ddgbs_for_rfid:
            filtered_rfids.append(rfid)
            selected_vcfs.append(last_ddgbs_for_rfid[rfid])
            ddgbs_selected += 1
        elif rfid in lsw_rnaseq_for_rfid:
            filtered_rfids.append(rfid)
            selected_vcfs.append(lsw_rnaseq_for_rfid[rfid])
            lsw_selected += 1
        else:
            # Exclude if no known pattern matches for this RFID.
            continue

    print(
        f"Parsed {len(vcf_list)} VCF IDs; Riptide RFIDs: {len(best_riptide_for_rfid)}, "
        f"ddGBS RFIDs: {len(last_ddgbs_for_rfid)}, LSWrnaseq RFIDs: {len(lsw_rnaseq_for_rfid)}; "
        f"unknown pattern lines: {unknown_patterns}; "
        f"selected {len(filtered_rfids)} of {len(req_list)} requested "
        f"(Riptide={riptide_selected}, ddGBS={ddgbs_selected}, LSWrnaseq={lsw_selected})",
        file=sys.stderr,
    )

    return filtered_rfids, selected_vcfs


def main(argv: Optional[List[str]] = None) -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Subset RFIDs to those with genotypes and select one VCF ID per RFID "
            "(choosing the last occurrence when multiple exist)."
        )
    )
    ap.add_argument("--rfids", required=True, help="Path to input RFIDs list (one per line)")
    ap.add_argument("--vcf-ids", required=True, help="Path to input VCF sample IDs list")
    ap.add_argument(
        "--out-tsv",
        required=True,
        help="Output path for TSV pairs: VCF_ID<TAB>RFID (no header)",
    )
    args = ap.parse_args(argv)

    with open(args.rfids, "r", encoding="utf-8") as f_in_r:
        rfids = f_in_r.read().splitlines()
    with open(args.vcf_ids, "r", encoding="utf-8") as f_in_v:
        vcf_ids = f_in_v.read().splitlines()

    filtered_rfids, selected_vcfs = select_vcf_per_rfid(rfids, vcf_ids)

    # Write one TSV with VCF_ID<TAB>RFID pairs (no header)
    with open(args.out_tsv, "w", encoding="utf-8") as f_out:
        for sid, rfid in zip(selected_vcfs, filtered_rfids):
            f_out.write(f"{sid}\t{rfid}\n")

    print(
        f"Wrote {len(selected_vcfs)} VCF/RFID pairs to {args.out_tsv}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
