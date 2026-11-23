"""Subset the HPA bulk TSVs to the biotin-ACC-FASN axis genes.

The bulk archives live in 02_Data/HPA (gitignored: re-downloadable from
proteinatlas.org, see 02_Data/HPA/urls.txt) and the extracted subsets are also
gitignored by the repo's 02_Data/**/*.tsv rule, so this script is the single
reproducible step between them and the figures.

Usage: python extract_hpa_axis.py   (needs the zips in 02_Data/HPA/tsv/)
"""
import os
import zipfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
HPA = os.path.join(HERE, os.pardir, os.pardir, "02_Data", "HPA")
ZIP_DIR = HPA          # the archives land here (aria2 -d), not in tsv/
TSV_DIR = os.path.join(HPA, "tsv")   # unzipped copies
OUT = os.path.join(HPA, "extract")

GENES = {
    "SLC5A6": "ENSG00000138074",   # SMVT itself; matches UniProt Q9Y289 / AF-Q9Y289-F1
    "ACACA": "ENSG00000278540",    # ACC1
    "ACACB": "ENSG00000076555",    # ACC2
    "FASN": "ENSG00000169710",
    "HLCS": "ENSG00000159267",     # holocarboxylase synthetase
    "BTD": "ENSG00000169814",      # biotinidase
}
FILES = ["normal_ihc_data.tsv", "cancer_data.tsv", "ms_tissue.tsv",
         "cancer_cptac.tsv", "subcellular_location.tsv", "proteinatlas.tsv"]


def main():
    os.makedirs(OUT, exist_ok=True)
    ids = set(GENES.values())
    for name in FILES:
        zpath = os.path.join(ZIP_DIR, name + ".zip")
        if not os.path.exists(zpath):
            print(f"SKIP {name}: {zpath} missing (download via 02_Data/HPA/urls.txt)")
            continue
        with zipfile.ZipFile(zpath) as z:
            with z.open(name) as fh:
                df = pd.read_csv(fh, sep="\t", dtype=str, low_memory=False)
        # proteinatlas.tsv keys genes by symbol, the others by Ensembl id
        col = "Gene" if df["Gene"].isin(ids).any() else "Ensembl"
        sub = df[df[col].isin(ids)] if col in df.columns else df[df["Gene"].isin(GENES)]
        out = os.path.join(OUT, name.replace(".tsv", "") + "_axis.tsv")
        sub.to_csv(out, sep="\t", index=False)
        print(f"{name:26s} total={len(df):7d} -> axis={len(sub):5d}  {out}")


if __name__ == "__main__":
    main()
