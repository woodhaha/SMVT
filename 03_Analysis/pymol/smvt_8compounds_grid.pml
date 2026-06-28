# ═══════════════════════════════════════════════════════
# SMVT MD Experiment — 8 Compounds Docking Grid
# All paths absolute for Windows PyMOL compatibility
# ═══════════════════════════════════════════════════════

# ── Clean slate ──
reinitialize
bg_color white
set valence, off
set grid_mode, on
set antialias, 2
set depth_cue, on

# ═══ Receptor ═══
load D:/Researching/SMVT/03_Analysis/pymol/SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey85, SMVT
set cartoon_transparency, 0.15, SMVT

# ── Pocket highlight ──
select pocket, SMVT and resid 116-142+228-265+309-331+350-378
show sticks, pocket
color lightorange, pocket
set stick_radius, 0.12, pocket
show surface, pocket
color grey95, pocket
set transparency, 0.55, pocket

# ═══ Load 8 ligands ═══
load D:/Researching/SMVT/03_Analysis/pymol/Hydromorphone_docked.pdb, Hydromorphone
load D:/Researching/SMVT/03_Analysis/pymol/Furosemide_docked.pdb, Furosemide
load D:/Researching/SMVT/03_Analysis/pymol/Naftazone_docked.pdb, Naftazone
load D:/Researching/SMVT/03_Analysis/pymol/Phenobarbital_docked.pdb, Phenobarbital
load D:/Researching/SMVT/03_Analysis/pymol/Biotin_docked.pdb, Biotin
load D:/Researching/SMVT/03_Analysis/pymol/Gabapentin_enacarbil_docked.pdb, Gabapentin
load D:/Researching/SMVT/03_Analysis/pymol/Riboflavin_docked.pdb, Riboflavin
load D:/Researching/SMVT/03_Analysis/pymol/Esketamine_docked.pdb, Esketamine

# ═══ Clone receptor 8x and translate ═══
create R1, SMVT
create R2, SMVT
create R3, SMVT
create R4, SMVT
create R5, SMVT
create R6, SMVT
create R7, SMVT
create R8, SMVT

# Row 1
translate [0, 0, 0], R1
translate [70, 0, 0], R2
translate [140, 0, 0], R3
translate [210, 0, 0], R4
# Row 2
translate [0, 80, 0], R5
translate [70, 80, 0], R6
translate [140, 80, 0], R7
translate [210, 80, 0], R8

# Translate ligands to match their receptors
translate [0, 0, 0], Hydromorphone
translate [70, 0, 0], Furosemide
translate [140, 0, 0], Naftazone
translate [210, 0, 0], Phenobarbital
translate [0, 80, 0], Biotin
translate [70, 80, 0], Gabapentin
translate [140, 80, 0], Riboflavin
translate [210, 80, 0], Esketamine

# ═══ Style ligands ═══
# Row 1 — Test Compounds
hide everything, Hydromorphone
show sticks, Hydromorphone
color red, Hydromorphone
set stick_radius, 0.25

hide everything, Furosemide
show sticks, Furosemide
color orange, Furosemide
set stick_radius, 0.25

hide everything, Naftazone
show sticks, Naftazone
color blue, Naftazone
set stick_radius, 0.25

hide everything, Phenobarbital
show sticks, Phenobarbital
color purple, Phenobarbital
set stick_radius, 0.25

# Row 2 — Controls + Pilot
hide everything, Biotin
show sticks, Biotin
color green, Biotin
set stick_radius, 0.25

hide everything, Gabapentin
show sticks, Gabapentin
color teal, Gabapentin
set stick_radius, 0.25

hide everything, Riboflavin
show sticks, Riboflavin
color grey70, Riboflavin
set stick_radius, 0.25

hide everything, Esketamine
show sticks, Esketamine
color cyan, Esketamine
set stick_radius, 0.25

# ═══ Labels (pseudoatom + label) ═══
set label_size, 18

pseudoatom L1, pos=[0, -28, 0]
label L1, "1. HYDROMORPHONE  -8.58 #1 Hit"
set label_color, red, L1

pseudoatom L2, pos=[70, -28, 0]
label L2, "2. FUROSEMIDE  -8.36 #2 Hit"
set label_color, orange, L2

pseudoatom L3, pos=[140, -28, 0]
label L3, "3. NAFTAZONE  -8.34 #3 Hit"
set label_color, blue, L3

pseudoatom L4, pos=[210, -28, 0]
label L4, "4. PHENOBARBITAL  -8.30 #4 Hit"
set label_color, purple, L4

pseudoatom L5, pos=[0, 108, 0]
label L5, "5. BIOTIN  -6.76  Natural Substrate"
set label_color, green, L5

pseudoatom L6, pos=[70, 108, 0]
label L6, "6. GABAPENTIN ENAC  -6.63  FDA Control"
set label_color, teal, L6

pseudoatom L7, pos=[140, 108, 0]
label L7, "7. RIBOFLAVIN  -0.01  Negative Ctrl"
set label_color, grey70, L7

pseudoatom L8, pos=[210, 108, 0]
label L8, "8. ESKETAMINE  -7.58  Pilot"
set label_color, cyan, L8

# Section headers
set label_size, 28
pseudoatom SEC1, pos=[105, 15, 0]
label SEC1, "<< TEST COMPOUNDS >>"
set label_color, grey30, SEC1

pseudoatom SEC2, pos=[105, 95, 0]
label SEC2, "<< CONTROLS + PILOT >>"
set label_color, grey30, SEC2

# ═══ Hide original SMVT ═══
hide everything, SMVT

# ═══ Zoom to see all ═══
zoom visible
extend visible, 1.3
move z, 60

# ═══ Print legend ═══
print "======================================================"
print " SMVT MD Experiment — 8 Compound Docking Grid"
print "======================================================"
print " TOP:    Hydromorphone | Furosemide | Naftazone | Phenobarbital"
print " BOTTOM: Biotin | Gabapentin | Riboflavin | Esketamine"
print "======================================================"
