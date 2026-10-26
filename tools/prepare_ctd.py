"""Extract bytewise sorted unique human chemical-gene pairs from the public CTD TSV."""
import argparse
import csv
import gzip
from pathlib import Path

def extract(source, target):
    seen=set()
    target.parent.mkdir(parents=True,exist_ok=True)
    with gzip.open(source,'rt',encoding='utf-8') as f:
        for row in csv.reader((line for line in f if not line.startswith('#')),delimiter='\t'):
            if len(row) < 11 or row[7] != '9606' or not row[0] or not row[3]:
                continue
            pair=(row[0],row[3])
            seen.add(pair)
    # Portable bytewise ordering; the analyzed null gene order is recorded separately.
    with target.open('w',encoding='utf-8',newline='') as out:
        writer=csv.writer(out,delimiter='\t',lineterminator='\n',quoting=csv.QUOTE_NONE)
        writer.writerows(sorted(seen,key=lambda p:(p[0].encode('utf-8'),p[1].encode('utf-8'))))
    return len(seen)

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    a=parser.parse_args()
    folder=a.root/'data'/'ctd'
    print('Human chemical-gene pairs:',extract(folder/'CTD_chem_gene_ixns.tsv.gz',folder/'ctd_human_chem_gene.tsv'))
