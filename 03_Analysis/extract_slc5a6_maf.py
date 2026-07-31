"""Extract SLC5A6 somatic mutations from full MC3 MAF into a small CSV (one-off)."""
import csv, sys

SRC = r'D:/Researching/SMVT/02_Data/raw/mc3.v0.2.8.PUBLIC.maf'
OUT = r'D:/Researching/SMVT/02_Data/cleaned/slc5a6_mc3_maf.csv'

rows = []
with open(SRC, 'r', encoding='utf-8') as f:
    reader = csv.reader(f, delimiter='\t')
    h = next(reader)
    for i, row in enumerate(reader):
        if row and row[0] == 'SLC5A6':
            rows.append(row)
        if i % 500000 == 0:
            sys.stderr.write(f'  scanned {i} lines, {len(rows)} SLC5A6 rows\n')

with open(OUT, 'w', encoding='utf-8', newline='') as f:
    w = csv.writer(f, delimiter='\t')
    w.writerow(h)
    w.writerows(rows)

sys.stderr.write(f'DONE: {len(rows)} SLC5A6 rows -> {OUT}\n')
