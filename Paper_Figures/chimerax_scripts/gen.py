#!/usr/bin/env python3
"""Generate ready-to-run ChimeraX .cxc scripts for all 8 SMVT compounds."""
import os

PDB_DIR = r"D:\Researching\SMVT\SMVT-MD-Analysis\chimerax"
OUT_DIR = r"D:\Researching\SMVT\Paper_Figures\chimerax_scripts"
os.makedirs(OUT_DIR, exist_ok=True)

C = ['HYDROMORPHONE','NAFTAZONE','PHENOBARBITAL','FUROSEMIDE',
     'BIOTIN','ESKETAMINE','GABAPENTIN_ENACARBIL','RIBOFLAVIN']
LAB = ['Hydromorphone','Naftazone','Phenobarbital','Furosemide',
       'Biotin','Esketamine','Gabapentin Enacarbil','Riboflavin']

for c, lbl in zip(C, LAB):
    pdb = os.path.join(PDB_DIR, f'{c}_centroid.pdb')
    out_png = os.path.join(OUT_DIR, f'{c}_binding_3d.png')
    script = f"""open {pdb}
hide solvent
cartoon; color lightgray target c
select :LIG; show sel; style sel stick; color sel byelement target a
~select sel
select zone :LIG range 5; show sel; style sel stick; color sel mediumpurple target a
~select sel
label /A name CA font helvetica 24 bold color black
findhbond selRestrict sel :LIG || selected reveal true dashes 12 color #22CC44 width 3
lighting full; set bgColor white; ~set glow
view orient
save {out_png} width 2400 height 1800 supersample 4 transparentBackground true
close session
"""
    with open(os.path.join(OUT_DIR, f'{c}.cxc'), 'w') as f:
        f.write(script)
    print(f'OK {c}.cxc')

# README
readme = f"""═══ SMVT ChimeraX 3D Binding Site Scripts ═══

How to use:
1. Open ChimeraX
2. File → Open → select any .cxc file
   OR drag-and-drop .cxc into ChimeraX window
3. Script auto-loads the centroid PDB, renders, saves PNG
4. PNG saved to same folder as .cxc

Compounds:
"""
for c, lbl in zip(C, LAB):
    readme += f"  {c}.cxc  →  {lbl}\n"
readme += f"""
Output: 2400×1800 PNG, transparent background, 4× supersample
PDBs from: {PDB_DIR}

Tip: To render all at once, open ChimeraX and run:
  cd {OUT_DIR}
  open *.cxc
"""
with open(os.path.join(OUT_DIR, 'README.txt'), 'w') as f:
    f.write(readme)
print('OK README.txt')
print(f'\nScripts in: {OUT_DIR}')
print('Open each .cxc in ChimeraX GUI to render 3D binding figures.')
