import csv
import argparse

def main():
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(description='Convert RefSeq accessions to chromosome names in GTF file')
    parser.add_argument('input_gtf', help='Input GTF file with RefSeq accessions')
    parser.add_argument('mapping_file', help='Tab-delimited file mapping RefSeq accessions to sequence names')
    parser.add_argument('output_gtf', help='Output GTF file with chromosome names')
    
    args = parser.parse_args()

    # Build a dictionary mapping from RefSeq accession to "chr" + Sequence name
    mapping = {}
    with open(args.mapping_file) as seq_file:
        reader = csv.DictReader(seq_file, delimiter="\t")
        for row in reader:
            refseq_acc = row["RefSeq seq accession"]
            seq_name = row["Sequence name"]
            if seq_name == "MT":
                seq_name = "M"
            mapping[refseq_acc] = "chr" + seq_name

    # Process the GTF file: update chromosome names in the first column
    with open(args.input_gtf) as gtf_in, open(args.output_gtf, "w") as gtf_out:
        for line in gtf_in:
            # Pass through header/comment lines unchanged
            if line.startswith("#"):
                gtf_out.write(line)
                continue
            fields = line.strip().split("\t")
            # Replace the chromosome ID if found in our mapping
            if fields[0] in mapping:
                fields[0] = mapping[fields[0]]
            gtf_out.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    main()
