"""Fetch Human Protein Atlas protein-level annotation for the SMVT lipogenic axis.

Why: the manuscript's SLC5A6 claims rest on bulk RNA-seq, plus single-cell in one
cancer type. HPA's per-gene record adds an orthogonal protein layer — whether
protein itself was detected, where it sits in the cell, and whether the protein
class flags it as a druggable transporter — and independently confirms the UniProt
accession the receptor model was built from.

Endpoint status verified on 2026-11-23 (older notes in this repo were stale):
  WORKS  https://www.proteinatlas.org/api/search_download.php?search=<GENE>
                 &format=json&columns=g,eg,p,up,pc,pe,rel,scl,scml,scal&compress=no
  DEAD   proteinatlas_external_dl.php        -> 404 "File not found."
  DEAD   subcellular_location_dl.php         -> 404 "File not found."
  DEAD   /<GENE>.tsv, format=tsv             -> 404 / 400
  NOTE   /ENSG<id>-<gene>-<any>.tsv serves the same ~269 kB HTML page whatever
         the data-type suffix is, so it is NOT a per-type TSV. Quantitative IHC
         and tissue nTPM therefore come from the bulk downloads under
         https://www.proteinatlas.org/download/tsv/ (normal_tissue_protein.tsv.gz,
         cancer_tissue_protein.tsv.gz, consensus_protein...), a separate step.

Output: 03_Analysis/outputs/hpa_protein/axis_annotation.json + axis_annotation.csv
"""
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request

GENES = ["SLC5A6", "ACACA", "ACACB", "FASN", "HLCS", "BTD", "SLC5A7"]
COLUMNS = "g,eg,p,up,pc,pe,rel,scl,scml,scal"
URL = "https://www.proteinatlas.org/api/search_download.php"
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "hpa_protein")


def fetch_gene(gene):
    qs = urllib.parse.urlencode({"search": gene, "format": "json",
                                 "columns": COLUMNS, "compress": "no"})
    req = urllib.request.Request(URL + "?" + qs,
                                 headers={"User-Agent": "Mozilla/5.0 (research)"})
    with urllib.request.urlopen(req, timeout=60) as r:
        body = r.read().decode("utf-8", "replace")
    time.sleep(1.2)  # HPA returns empty responses above ~4 req/s
    return json.loads(body)


def main():
    os.makedirs(OUT, exist_ok=True)
    records = {}
    for gene in GENES:
        rows = fetch_gene(gene)
        exact = [r for r in rows if r.get("Gene") == gene]
        if not exact:
            print(f"WARN {gene}: no exact match among {len(rows)} rows")
            continue
        rec = exact[0]
        records[gene] = rec
        print(f"{gene:8s} {rec.get('Ensembl','?'):18s} UniProt={','.join(rec.get('Uniprot') or [])[:10]:10s} "
              f"{rec.get('Evidence','?')} | {','.join(rec.get('Subcellular main location') or [])}")
    with open(os.path.join(OUT, "axis_annotation.json"), "w", encoding="utf-8") as fh:
        json.dump(records, fh, indent=1, ensure_ascii=False)
    flat_fields = ["Gene", "Ensembl", "Uniprot", "Evidence", "Protein class",
                   "Reliability score", "Subcellular location",
                   "Subcellular main location", "Subcellular additional location"]
    with open(os.path.join(OUT, "axis_annotation.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Field"] + list(records))
        for field in flat_fields:
            values = []
            for g in records:
                v = records[g].get(field)
                values.append("; ".join(v) if isinstance(v, list) else (v or ""))
            w.writerow([field] + values)
    print(f"\nwrote {len(records)} genes -> {OUT}")
    return 0 if records else 1


if __name__ == "__main__":
    sys.exit(main())
