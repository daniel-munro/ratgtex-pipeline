"""Compare the actual chromosome-filtered RatGTEx reference annotations.
Run from the RatGTEx repository root: python reports/annotation_v4_vs_v5/compare.py
"""
import collections
import csv
import json
import pathlib
import re

OUT = pathlib.Path(__file__).resolve().parent
PATHS = {
    'v4': 'v4/ref/GCF_036323735.1_GRCr8_genomic.chr.gtf',
    'v5': 'ref/GCF_036323735.1-RS_2026_08_genomic.chr.gtf',
}

def load(path):
    genes, tx, exons, cds, headers = {}, {}, set(), set(), []
    gene_records = 0
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                headers.append(line.strip())
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = line.rstrip().split('\t')
            if feature not in {'gene', 'transcript', 'exon', 'CDS'}:
                continue
            a = dict(re.findall(r'(\w+) "([^"]*)"', attrs))
            gid = re.search(r'GeneID:(\d+)', attrs).group(1)
            if feature == 'gene':
                gene_records += 1
                record = dict(symbol=a['gene_id'], biotype=a['gene_biotype'], chrom=chrom, start=int(start), end=int(end), strand=strand)
                if gid in genes:
                    assert genes[gid]['biotype'] == record['biotype'], (gid, path)
                    genes[gid]['locations'].add((chrom, int(start), int(end), strand))
                else:
                    genes[gid] = dict(**record, locations={(chrom, int(start), int(end), strand)})
            elif feature == 'transcript':
                tid = a['transcript_id']
                assert tid not in tx or tx[tid] == gid, (tid, path)
                tx[tid] = gid
            elif feature == 'exon':
                exons.add((chrom, start, end, strand))
            elif feature == 'CDS':
                cds.add(gid)
    counts = collections.Counter(g['biotype'] for g in genes.values())
    tx_per_gene = collections.Counter(tx.values())
    stats = dict(gene_feature_records=gene_records, genes=len(genes), transcripts=len(tx), unique_exon_intervals=len(exons), genes_with_CDS=len(cds), genes_with_transcripts=len(tx_per_gene), genes_with_multiple_transcripts=sum(n > 1 for n in tx_per_gene.values()), **{'genes:'+k:v for k,v in counts.items()})
    return genes, tx, stats, headers

data = {v:load(p) for v,p in PATHS.items()}
with open(OUT/'counts.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['metric','v4','v5','difference','percent_change'])
    for key in sorted(set(data['v4'][2]) | set(data['v5'][2])):
        a,b = [data[v][2].get(key,0) for v in PATHS]
        w.writerow([key,a,b,b-a,round(100*(b-a)/a,2) if a else 'NA'])
g4,g5 = data['v4'][0],data['v5'][0]
shared = g4.keys() & g5.keys()
comparison = dict(shared_GeneIDs=len(shared), v4_only_GeneIDs=len(g4.keys()-g5.keys()), v5_only_GeneIDs=len(g5.keys()-g4.keys()), shared_genes_renamed=sum(g4[g]['symbol']!=g5[g]['symbol'] for g in shared), shared_genes_changed_biotype=sum(g4[g]['biotype']!=g5[g]['biotype'] for g in shared), shared_genes_changed_location=sum(g4[g]['locations']!=g5[g]['locations'] for g in shared))
with open(OUT/'gene_changes.tsv','w') as f:
    w=csv.writer(f,delimiter='\t'); w.writerow(['GeneID','status','v4_symbol','v5_symbol','v4_biotype','v5_biotype'])
    for g in sorted(g4.keys()|g5.keys(),key=int):
        a,b=g4.get(g),g5.get(g)
        w.writerow([g,'shared' if a and b else 'v4_only' if a else 'v5_only',a['symbol'] if a else '',b['symbol'] if b else '',a['biotype'] if a else '',b['biotype'] if b else ''])
with open(OUT/'provenance.json','w') as f:
    json.dump(dict(inputs={v:dict(path=p,headers=data[v][3]) for v,p in PATHS.items()},comparison=comparison),f,indent=2)
print(json.dumps(comparison,indent=2))
print((OUT/'counts.tsv').read_text())
