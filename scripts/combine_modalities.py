"""Combine all configured modalities, requiring identical sample headers."""
import gzip
import subprocess
import sys
from pathlib import Path

project = Path(sys.argv[1])
modalities = sys.argv[2:]
output = project / 'phenos/output'
header = None
with (output / 'cross_modality.unsorted.bed').open('w') as bed, (output / 'cross_modality.phenotype_groups.txt').open('w') as groups:
    for modality in modalities:
        with gzip.open(output / f'{modality}.bed.gz', 'rt') as handle:
            current_header = next(handle)
            if header is None:
                header = current_header
            elif header != current_header:
                raise ValueError(f'Sample header differs for {modality}')
            for line in handle:
                fields = line.rstrip('\n').split('\t')
                fields[3] = f'{modality}:{fields[3]}'
                bed.write('\t'.join(fields) + '\n')
                groups.write(f"{fields[3]}\t{fields[3].split(':', 1)[1].split('__', 1)[0]}\n")
with (output / 'cross_modality.bed').open('w') as handle:
    handle.write(header)
    handle.flush()  # Write the header before sort writes directly to this file descriptor.
    subprocess.run(['sort', '-T', str(output), '-k1,1', '-k2,2n', str(output / 'cross_modality.unsorted.bed')], stdout=handle, check=True)
subprocess.run(['bgzip', '-f', str(output / 'cross_modality.bed')], check=True)
subprocess.run(['tabix', '-f', '-p', 'bed', str(output / 'cross_modality.bed.gz')], check=True)
(output / 'cross_modality.unsorted.bed').unlink()
