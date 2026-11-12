"""Render designed molecules with functional-group highlighting (RDKit only, no matplotlib
to avoid a Cairo/DLL clash in the smvt-md env). Legend is drawn by make_legend.py."""
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

FIG = r"D:\Researching\SMVT\autoresearch\paper\figures"
RED = (0.714, 0.263, 0.259)    # COOH
TEAL = (0.259, 0.580, 0.620)   # CF3
GREEN = (0.180, 0.620, 0.267)  # linker


def smarts_atoms(m, sma):
    p = Chem.MolFromSmarts(sma)
    out = set()
    for match in m.GetSubstructMatches(p):
        out.update(match)
    return sorted(out)


MOLS = [
    ("COOH-bisCF3", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2",
     "COOH-bisCF3 (parent)", None),
    ("COOH-bisCF3-NH2", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CN)",
     "COOH-bisCF3-NH2 (amine handle)", "[NX3;H2]"),
    ("COOH-bisCF3-alkyne", "OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CC#C)",
     "COOH-bisCF3-alkyne (click handle)", "C#C"),
    ("COOH-bisCF3-SN38", "CCC1=C2C=C(OC(=O)NCC34CC5(C(F)(F)F)CC(CC(C5)(C3)C(F)(F)F)C4C(=O)O)C=CC2=NC2=C1CN1C(=O)C3=C(C=C21)C(CC)(O)COC(=O)3",
     "COOH-bisCF3-SN38 (carbamate prodrug)", "[NX3;H0][CX3](=O)[OX2]"),
]

for name, smi, desc, linker_sma in MOLS:
    m = Chem.MolFromSmiles(smi)
    AllChem.Compute2DCoords(m)
    cooh = smarts_atoms(m, "[CX3](=O)[OX2H1]") + smarts_atoms(m, "[OX1]=[CX3]")
    cf3 = smarts_atoms(m, "C(F)(F)F")
    linker = smarts_atoms(m, linker_sma) if linker_sma else []
    highlight = cooh + cf3 + linker
    colors = {}
    for a in cooh:
        colors[a] = RED
    for a in cf3:
        colors[a] = TEAL
    for a in linker:
        colors[a] = GREEN
    img = Draw.MolToImage(m, size=(720, 480), highlightAtoms=highlight,
                          highlightAtomColors=colors,
                          highlightAtomRadii={a: 0.35 for a in highlight})
    img.save(FIG + "/" + name + "_annotated.png")
    print("saved", name)
