#!/usr/bin/env python3

"""Summarize STAR Log.final.out statistics for multiple samples."""

import argparse
from pathlib import Path

import pandas as pd

def parse_star_log(log_file):
    """Parse a STAR Log.final.out file into a dict of label -> value."""
    stats = {}
    with open(log_file, 'r') as f:
        for line in f:
            line = line.strip()
            if '|' in line:
                label, value = line.split('|')
                label = label.strip()
                value = value.strip()
                stats[label] = value
    return stats

def main():
    parser = argparse.ArgumentParser(description='Extract statistics from STAR log files')
    parser.add_argument('star_dir', type=Path, help='Directory containing STAR log files')
    parser.add_argument('samples_file', type=Path, help='File containing sample IDs')
    parser.add_argument('output_tsv', type=Path, help='Output TSV file path')
    
    args = parser.parse_args()
    
    with open(args.samples_file, 'r') as f:
        sample_ids = [line.strip() for line in f]
    
    # Parse each log file
    all_stats = []
    for sample_id in sample_ids:
        log_file = args.star_dir / f"{sample_id}.Log.final.out"
            
        stats = parse_star_log(log_file)
        stats = {'sample_id': sample_id, **stats}  # Put Sample first by using dict unpacking
        all_stats.append(stats)
    
    df = pd.DataFrame(all_stats)
    df.to_csv(args.output_tsv, sep='\t', index=False)

if __name__ == '__main__':
    main() 
