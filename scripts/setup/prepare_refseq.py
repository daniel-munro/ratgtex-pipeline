import csv
import argparse


def load_mapping(mapping_file: str) -> dict:
    """Load mapping from RefSeq accession to common chromosome name.

    The mapping file is expected to be a TSV with headers including:
      - "RefSeq seq accession"
      - "Sequence name"
    """
    mapping: dict[str, str] = {}
    with open(mapping_file) as seq_file:
        reader = csv.DictReader(seq_file, delimiter="\t")
        for row in reader:
            refseq_acc = row.get("RefSeq seq accession")
            seq_name = row.get("Sequence name")
            mapping[refseq_acc] = seq_name
    return mapping


def process_gtf(input_gtf: str, output_gtf: str, mapping: dict) -> None:
    """Update the first column (seqname) of a GTF using the provided mapping."""
    with open(input_gtf) as gtf_in, open(output_gtf, "w") as gtf_out:
        for line in gtf_in:
            # Pass through header/comment lines unchanged
            if line.startswith("#"):
                gtf_out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            # Replace the chromosome ID if found in our mapping
            if fields[0] in mapping:
                fields[0] = mapping[fields[0]]
            gtf_out.write("\t".join(fields) + "\n")


def process_fasta(input_fasta: str, output_fasta: str, mapping: dict) -> None:
    """Update FASTA headers by replacing the first token (RefSeq accession)
    with the mapped chr name; preserve the rest of the header and all sequences.
    """
    with open(input_fasta) as fin, open(output_fasta, "w") as fout:
        for line in fin:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                # Split once on the first space to separate ID and the rest
                if " " in header:
                    acc, rest = header.split(" ", 1)
                else:
                    acc, rest = header, ""
                new_id = mapping.get(acc, acc)
                if rest:
                    fout.write(f">{new_id} {rest}\n")
                else:
                    fout.write(f">{new_id}\n")
            else:
                fout.write(line)

def main():
    parser = argparse.ArgumentParser(
        description="Convert RefSeq accessions to chromosome names in GTF and optionally FASTA files"
    )
    parser.add_argument('--map', help='Tab-delimited file mapping RefSeq accessions to sequence names')
    parser.add_argument('--gtf-in', help='Input GTF file with RefSeq accessions')
    parser.add_argument('--gtf-out', help='Output GTF file with chromosome names')
    parser.add_argument('--fasta-in', help='Input reference genome FASTA with RefSeq accessions in headers')
    parser.add_argument('--fasta-out', help='Output FASTA with chromosome names in headers')
    
    args = parser.parse_args()

    mapping = load_mapping(args.map)
    process_gtf(args.gtf_in, args.gtf_out, mapping)
    if args.fasta_in and args.fasta_out:
        process_fasta(args.fasta_in, args.fasta_out, mapping)

if __name__ == "__main__":
    main()
