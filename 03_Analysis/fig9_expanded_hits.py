"""Fig. 9 | Expanded-library screening hits (20 compounds with dG < -8.0 kcal/mol).

Rebuilds a figure that previously existed only as an oversized PDF (449 x 254 mm,
2.5x the 183 mm journal maximum, no producing script): drawn at double-column width
with labels above the 5 pt print floor, and emitted through fig_style.save_pub so the
PDF/SVG/PNG/LZW-TIFF set is reproducible from the docking output alone.

Data: 03_Analysis/outputs/docking_final/hits20_for_fig.json (name, dG, chemotype,
smiles). Structures are RDKit 2D depictions of those SMILES; dG values are printed
verbatim from the file, so the figure cannot drift from the docking results.
"""
import io
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
from PIL import Image
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

sys.path.insert(0, r"D:\Researching\SMVT\03_Analysis")
from fig_style import apply_style, save_pub, DOUBLE, check_text_floor

SRC = r"D:\Researching\SMVT\03_Analysis\outputs\docking_final\hits20_for_fig.json"
OUTDIR = r"D:\Researching\SMVT\04_Manuscript\figures"
NCOL, NROW = 5, 4
apply_style()

hits = json.load(open(SRC, encoding="utf-8"))
hits = sorted(hits, key=lambda h: float(h["dG"]))          # strongest first, left-to-right
print(f"{len(hits)} compounds, dG from {hits[0]['dG']} to {hits[-1]['dG']}")

fig, axes = plt.subplots(NROW, NCOL, figsize=(DOUBLE, 150 / 25.4))
for ax, h in zip(axes.ravel(), hits):
    ax.axis("off")
    mol = Chem.MolFromSmiles(h["smiles"])
    if mol is None:
        ax.text(0.5, 0.5, f"unparsable\n{h['name'][:18]}", ha="center", va="center", fontsize=5)
        continue
    d = rdMolDraw2D.MolDraw2DCairo(360, 300)
    d.drawOptions().addStereoAnnotation = False
    d.DrawMolecule(mol)
    d.FinishDrawing()
    img = Image.open(io.BytesIO(d.GetDrawingText()))
    ax.imshow(img, interpolation="bilinear")
    ax.set_title(f"{h['name'][:24]}{'..' if len(h['name']) > 24 else ''}\n"
                 f"{float(h['dG']):.2f} kcal/mol   {h['chemotype']}",
                 fontsize=5.6, pad=2)

fig.suptitle("Fig. 9 | Expanded-library hits (dG < -8.0 kcal/mol), strongest first",
             fontsize=7.5, fontweight="bold", y=0.995)
plt.tight_layout(rect=(0, 0, 1, 0.975))

low = check_text_floor(fig, print_mm=183)
files = save_pub(fig, os.path.join(OUTDIR, "Fig9_expanded_screening_hits"))
print("wrote:", [os.path.basename(f) for f in files])
print("sub-5pt text:", low if low else "none")
print("\ndG values now on the figure (check these against the manuscript text):")
for h in hits:
    print(f"  {float(h['dG']):+.2f}  {h['name'][:44]:44s} {h['chemotype']}")
plt.close(fig)
