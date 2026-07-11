"""
PyMOL visualization script for SMVT binding modes.
Usage: pymol -cq _pymol_vis.py
"""
import sys, os

# Configuration
ANALYSIS = r'D:\Researching\SMVT\MM-GBSA Analysis'
TOP_COMPOUNDS = ['HYDROMORPHONE', 'BIOTIN', 'GABAPENTIN_ENACARBIL', 'NAFTAZONE']

# Colors (PyMOL names)
COLORS = {
    'BIOTIN': 'marine',
    'HYDROMORPHONE': 'purple',
    'GABAPENTIN_ENACARBIL': 'lime',
    'NAFTAZONE': 'orange',
}

# Start PyMOL
import pymol
pymol.finish_launching(['pymol', '-cq'])

# Load receptor (use one receptor, e.g. BIOTIN receptor)
cmd = pymol.cmd
cmd.load(os.path.join(ANALYSIS, 'BIOTIN_receptor.pdb'), 'receptor')

# Show receptor as cartoon
cmd.hide('everything', 'receptor')
cmd.show('cartoon', 'receptor')
cmd.color('white', 'receptor')

# Select and show binding pocket (conserved residues from pocket_contacts.json)
# Based on consensus residues found in all 8 compounds
pocket_res = 'resid 76+79+80+81+84+88+99+106+156+259+260+263+264+265+266+267+268+270+271+277+297+300+301+302+305+366+424+425+427+428+431+432+527+528'
cmd.select('pocket', f'receptor and chain A and {pocket_res}')
cmd.show('sticks', 'pocket')
cmd.color('palegreen', 'pocket')
cmd.select('pocket_ca', f'receptor and chain A and {pocket_res} and name CA')
cmd.show('spheres', 'pocket_ca')
cmd.color('green', 'pocket_ca')

# Highlight key H-bonding residues
cmd.select('hbond_res', 'receptor and chain A and resi 79+80+84+99+106+271+277+301+425+431')
cmd.show('sticks', 'hbond_res')
cmd.color('cyan', 'hbond_res')
cmd.label('hbond_res and name CA', 'resn + " " + resi')

# Load docked poses of top compounds
for name in TOP_COMPOUNDS:
    fn = os.path.join(ANALYSIS, f'{name}_docked.pdb')
    if not os.path.exists(fn):
        print(f'{name}: PDB not found at {fn}')
        continue
    cmd.load(fn, name)
    cmd.hide('everything', name)
    cmd.show('sticks', name)
    cmd.color(COLORS.get(name, 'gray'), name)
    cmd.label(f'{name} and name C1', '"' + name + '"')
    cmd.set('label_size', 1.0)

# Set view
cmd.orient('pocket')
cmd.zoom('pocket', 5)

# Create figure
img_path = os.path.join(ANALYSIS, 'ALL_pymol_binding.png')
cmd.png(img_path, width=1600, height=1200, dpi=150, ray=1)
print(f'Saved: {img_path}')

# Also save session
sess_path = os.path.join(ANALYSIS, 'SMVT_binding.pse')
cmd.save(sess_path)
print(f'Saved: {sess_path}')

cmd.quit()
print('DONE')
