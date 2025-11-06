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
TISSUES = config["run"]
# Get list of tissues pre-merge for things like QC
TISSUES_SEP = [t for t in TISSUES if t not in config["merged_tissues"]]
TISSUES_MERGED = [t for t in TISSUES if not t[-1].isdigit()]
FASTQ_MAPS = {}
for tissue in TISSUES:
    if tissue in config["merged_tissues"]:
        config["tissues"][tissue] = {}
        dsets = config["merged_tissues"][tissue]
        geno_datasets = [config["tissues"][dset]["geno_dataset"] for dset in dsets]
        assert len(set(geno_datasets)) == 1, f"All merged tissues must have the same geno_dataset: {tissue}"
        config["tissues"][tissue]["geno_dataset"] = geno_datasets[0]
    else:
        fastq_path = Path(config["tissues"][tissue]["fastq_path"])
        paired_end = config["tissues"][tissue].get("paired_end", None)
        FASTQ_MAPS[tissue] = load_fastq_map(f"{VERSION}/{tissue}/fastq_map.txt", fastq_path, paired_end)
