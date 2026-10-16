"""Choose unique biallelic SNP positions and orient alleles to the reference FASTA."""
from collections import Counter
import gzip
from itertools import groupby
import json

import pysam


def pvar_rows(path):
    with open(path) as handle:
        for line in handle:
            if line.startswith('##'):
                continue
            header = line.rstrip().split('\t')
            break
        if header[:5] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT']:
            raise ValueError(f'Unexpected PVAR header in {path}')
        for line in handle:
            fields = line.rstrip().split('\t')
            yield fields[:5]


def prepare_variants(pvar, fasta, chromosome, extract, reference, audit, summary):
    with pysam.FastaFile(str(fasta)) as genome:
        sequence = genome.fetch(f'chr{chromosome}').upper()
    counts = Counter()
    seen_ids = set()
    previous = 0
    with open(extract, 'w') as keep, open(reference, 'w') as refs, gzip.open(audit, 'wt') as report:
        report.write('chromosome\tposition\tinput_id\tref\talt\tfasta_ref\taction\n')
        for pos, group in groupby(pvar_rows(pvar), key=lambda row: int(row[1])):
            rows = list(group)
            if pos <= previous or pos > len(sequence):
                raise ValueError(f'Unsorted or invalid position {pos} in {pvar}')
            previous = pos
            base = sequence[pos - 1]
            for chrom, _, variant, ref, alt in rows:
                if chrom.removeprefix('chr') != str(chromosome):
                    raise ValueError(f'Unexpected chromosome {chrom} in {pvar}')
                if variant in seen_ids or variant == '.':
                    raise ValueError(f'Duplicate or missing input variant ID: {variant}')
                seen_ids.add(variant)
                if len(rows) > 1:
                    action = 'exclude_duplicate_position'
                elif ref not in 'ACGT' or alt not in 'ACGT' or len(ref) != 1 or len(alt) != 1 or ref == alt:
                    action = 'exclude_non_biallelic_snp'
                elif base == ref:
                    action = 'keep'
                elif base == alt:
                    action = 'swap_ref_alt'
                else:
                    action = 'exclude_reference_mismatch'
                counts[action] += 1
                if action in ('keep', 'swap_ref_alt'):
                    keep.write(variant + '\n')
                    refs.write(f'{variant}\t{base}\n')
                if action != 'keep':
                    report.write(f'{chrom}\t{pos}\t{variant}\t{ref}\t{alt}\t{base}\t{action}\n')
    with open(summary, 'w') as handle:
        json.dump(dict(chromosome=str(chromosome), input_pvar=str(pvar),
                       reference_fasta=str(fasta), counts=dict(counts)), handle, indent=2)
        handle.write('\n')
    if not counts['keep'] + counts['swap_ref_alt']:
        raise ValueError(f'No usable SNPs in {pvar}')


if __name__ == '__main__':
    prepare_variants(snakemake.input.pvar, snakemake.input.fasta, snakemake.wildcards.chrom,
                     snakemake.output.extract, snakemake.output.reference,
                     snakemake.output.audit, snakemake.output.summary)
