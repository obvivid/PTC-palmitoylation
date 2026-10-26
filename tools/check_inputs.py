"""Verify public input snapshots against the manuscript input fingerprints."""
import argparse
import csv
import hashlib
from pathlib import Path

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):
            h.update(block)
    return h.hexdigest()

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    a=parser.parse_args()
    manifest=Path(__file__).resolve().parents[1]/'config'/'input_sha256.tsv'
    contract=Path(__file__).resolve().parents[1]/'config'/'inputs.csv'
    optional=set()
    if contract.is_file():
        with contract.open(encoding='utf-8-sig',newline='') as f:
            optional={r['file'] for r in csv.DictReader(f) if not r['steps'].strip()}
    missing=mismatch=optional_missing=0
    with manifest.open(encoding='utf-8') as f:
        for row in csv.DictReader(f,delimiter='\t'):
            p=a.root/row['file']
            if not p.is_file() and row['file'] in optional:
                print('OPTIONAL MISSING',row['file']); optional_missing+=1
            elif not p.is_file():
                print('MISSING',row['file']); missing+=1
            elif digest(p)!=row['sha256']:
                print('DIFFERENT SNAPSHOT',row['file']); mismatch+=1
            else:
                print('MATCH',row['file'])
    print(f'Input check: {missing} required missing; {mismatch} different snapshots; {optional_missing} optional missing')
    raise SystemExit(1 if missing or mismatch else 0)
