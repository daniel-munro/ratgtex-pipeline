def load_fastq_map(map_file: Path, fastq_dir: Path, paired_end=None) -> dict:
    """Load the FASTQ file paths for each sample.
    
    Returns a dictionary where each value is a tuple (paths, is_paired).
    paths is either a list (single-end) or tuple of two lists (paired-end).
    is_paired is a boolean indicating if the sample is paired-end.
    """
    paths = {}
    with open(map_file, "r") as f:
        for line in f.read().splitlines():
            fields = line.split('\t')
            rat_id = fields[-1]
            is_paired = len(fields) == 3
            
            if rat_id not in paths:
                paths[rat_id] = (([], []), is_paired) if is_paired else ([], is_paired)
            elif paths[rat_id][1] != is_paired:
                raise ValueError(f"Sample {rat_id} has mixed paired-end and single-end reads")
            
            if is_paired:
                paths[rat_id][0][0].append(str(fastq_dir / fields[0]))
                paths[rat_id][0][1].append(str(fastq_dir / fields[1]))
            else:
                paths[rat_id][0].append(str(fastq_dir / fields[0]))
    
    # Validate against config if specified
    if paired_end is not None:
        for rat_id, (_, is_paired) in paths.items():
            if is_paired != paired_end:
                raise ValueError(f"Sample {rat_id} has {'paired' if is_paired else 'single'}-end reads, but config specifies {'paired' if paired_end else 'single'}-end")
    
    return paths


VERSION = config["version"]
GENOME_PREFIX = config["ref_genome"].removesuffix(".fa")
ANNO_PREFIX = config["ref_anno"].removesuffix(".gtf")
import sys
sys.path.insert(0, str(Path(workflow.basedir) / 'scripts'))
from workflow_config import datasets

DATASETS = datasets(config)
TISSUES = list(config['tissues'])
MODALITIES = list(config['modalities'])
grouped_modality = {m: s['grouped'] for m, s in config['modalities'].items()}
if config['cross_modality']:
    grouped_modality['cross_modality'] = True
FASTQ_MAPS = {}
for dataset in DATASETS:
    settings = config['datasets'][dataset]
    FASTQ_MAPS[dataset] = load_fastq_map(
        f"{VERSION}/{dataset}/fastq_map.txt", Path(settings['fastq_path']), settings.get('paired_end'))
GENOTYPE_PANELS = {d: settings['genotype_panel'] for d, settings in config['datasets'].items()}
for tissue, members in config['tissues'].items():
    genotype_panels = {config['datasets'][d]['genotype_panel'] for d in members}
    if len(genotype_panels) != 1:
        raise ValueError(f'All datasets in tissue {tissue} must share a genotype panel')
    GENOTYPE_PANELS[tissue] = genotype_panels.pop()
