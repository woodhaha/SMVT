"""
Convert top-hit docked PDBQT → PDB for PyMOL, then generate .pml script.
"""
import os, glob

os.chdir(os.path.join(os.path.dirname(__file__), ".."))
os.makedirs("03_Analysis/pymol", exist_ok=True)

TOP_HITS = {
    # Has docked files
    "Furosemide":    "#E74C3C",   # −8.36  sulfonamide, WHO essential
    "Naftazone":     "#2980B9",   # −8.34  naphthoquinone
    "Phenobarbital": "#8E44AD",   # −8.30  barbiturate (100% class hit)
    "Diclofenac":    "#E67E22",   # −7.15  NSAID (NSAID-SMVT axis)
    "Esketamine":    "#3498DB",   # −7.58  arylcyclohexylamine, novel
    "Biotin":        "#27AE60",   # −6.76  natural substrate (reference)
    # Hydromorphone (−8.58 best) docked separately — not in this batch
}

def pdbqt_to_pdb(pdbqt_path, out_path):
    """Extract MODEL 1 from PDBQT and write as HETATM PDB."""
    atoms = []
    in_model = False
    with open(pdbqt_path) as f:
        for line in f:
            if line.startswith("MODEL") and int(line.split()[1]) == 1:
                in_model = True
            elif line.startswith("ENDMDL") and in_model:
                break
            elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                atoms.append(line[:66].ljust(80))
    if atoms:
        with open(out_path, "w") as f:
            f.write("\n".join(atoms) + "\nEND\n")
    return len(atoms)

# Convert docked poses
for name in TOP_HITS:
    pdbqt = f"03_Analysis/docking/{name}_docked.pdbqt"
    out = f"03_Analysis/pymol/{name}_docked.pdb"
    if os.path.exists(pdbqt):
        n = pdbqt_to_pdb(pdbqt, out)
        print(f"  {name}: {n} atoms → {out}")

# Copy receptor
import shutil
receptor_src = "02_Data/cleaned/SMVT_prepared.pdb"
receptor_dst = "03_Analysis/pymol/SMVT_receptor.pdb"
if os.path.exists(receptor_src):
    shutil.copy(receptor_src, receptor_dst)
    print(f"  Receptor: {receptor_src} → {receptor_dst}")
elif os.path.exists("02_Data/raw/AF-Q9Y289-F1.pdb"):
    shutil.copy("02_Data/raw/AF-Q9Y289-F1.pdb", receptor_dst)
    print(f"  Receptor (raw): → {receptor_dst}")

# ═══ Write PyMOL script ═══
pml = []
pml.append("""# SMVT (SLC5A6) Docking Results — PyMOL Session
# Top 4 Hits + Biotin Control
# Usage: pymol 03_Analysis/pymol/smvt_docking.pml

# ── Setup ──
bg_color white
set valence, off
set cartoon_fancy_helices, on
set cartoon_highlight_color, grey50
set ray_trace_mode, 1
set ray_trace_gain, 0.1

# ── Load receptor ──
load 03_Analysis/pymol/SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey80, SMVT
set cartoon_transparency, 0.2, SMVT

# ── Binding pocket highlight ──
# Select residues within 5Å of any docked ligand (approximate)
select pocket, SMVT and resi 1-600
show surface, pocket
color grey90, pocket
set transparency, 0.7, pocket

""")

for name, color in TOP_HITS.items():
    pml.append(f"""
# ── {name} ──
load 03_Analysis/pymol/{name}_docked.pdb, {name}
hide everything, {name}
show sticks, {name}
color {color}, {name}
set stick_radius, 0.2, {name}
""")

# Label key features
pml.append("""
# ── View ──
zoom SMVT
turn x, 30
turn y, -15

# ── Key interactions (manual annotation) ──
# Biotin binding site: ureido ring in transmembrane cavity
# Barbiturates: malonylurea overlaps with biotin ureido
# Show H-bonds (if donor/acceptor within 3.5Å):
# distance hbonds, SMVT, Biotin, 3.5, mode=2

# ── Group objects ──
group Hits, Furosemide Naftazone Phenobarbital Diclofenac Esketamine
group Controls, Biotin

# ── Legends ──
set label_size, 20
set label_color, black

# ── Save session ──
# save 03_Analysis/pymol/smvt_docking.pse

# ── Render views ──
set ray_opaque_background, off
viewport 1200, 800

# Front view
zoom SMVT
ray 1200, 800
# png 03_Analysis/pymol/smvt_docking_front.png

# Side view
turn y, 90
ray 1200, 800
# png 03_Analysis/pymol/smvt_docking_side.png

# Top-down (binding pocket)
orient pocket
zoom pocket, 2
ray 1200, 800
# png 03_Analysis/pymol/smvt_docking_pocket.png

print("Session ready. Use: ray; png output.png")
print("Controls are in group 'Controls', hits in group 'Hits'")
""")

pml_path = "03_Analysis/pymol/smvt_docking.pml"
with open(pml_path, "w") as f:
    f.write("\n".join(pml))

print(f"\n✓ PyMOL script: {pml_path}")
print("  Run: pymol 03_Analysis/pymol/smvt_docking.pml")
