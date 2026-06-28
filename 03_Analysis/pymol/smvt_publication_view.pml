# ═══════════════════════════════════════════════════════
# SMVT — Publication-Ready Docking Visualization
# 4-panel: Best Hit | Substrate | Pocket Detail | Surface
# ═══════════════════════════════════════════════════════

reinitialize
bg_color white
set valence, off
set antialias, 3
set depth_cue, on
set ray_trace_mode, 1
set ray_trace_gain, 0
set cartoon_fancy_helices, on
set cartoon_highlight_color, grey40
set stick_quality, 15
set sphere_quality, 2
set surface_quality, 2

# ═══ Load all assets ═══
load D:/Researching/SMVT/03_Analysis/pymol/SMVT_receptor.pdb, SMVT
load D:/Researching/SMVT/03_Analysis/pymol/Hydromorphone_docked.pdb, HM
load D:/Researching/SMVT/03_Analysis/pymol/Biotin_docked.pdb, Biotin
load D:/Researching/SMVT/03_Analysis/pymol/Furosemide_docked.pdb, FUR
load D:/Researching/SMVT/03_Analysis/pymol/Naftazone_docked.pdb, NAF
load D:/Researching/SMVT/03_Analysis/pymol/Phenobarbital_docked.pdb, PHB
load D:/Researching/SMVT/03_Analysis/pymol/Gabapentin_enacarbil_docked.pdb, GBP
load D:/Researching/SMVT/03_Analysis/pymol/Riboflavin_docked.pdb, RIB
load D:/Researching/SMVT/03_Analysis/pymol/Esketamine_docked.pdb, ESK
load D:/Researching/SMVT/03_Analysis/pymol/Pentobarbital_docked.pdb, PTB
load D:/Researching/SMVT/03_Analysis/pymol/Carprofen_docked.pdb, CAR
load D:/Researching/SMVT/03_Analysis/pymol/Diclofenac_docked.pdb, DIC
load D:/Researching/SMVT/03_Analysis/pymol/Butabarbital_docked.pdb, BTB

# ═══ PANEL A: Hydromorphone (−8.58) — Best Hit ═══
create pA_receptor, SMVT
create pA_ligand, HM
translate [0, 0, 0], pA_receptor
translate [0, 0, 0], pA_ligand

# Receptor
hide everything, pA_receptor
show cartoon, pA_receptor
color grey85, pA_receptor

# Pocket residues (within 5A)
select pA_pocket, pA_receptor within 5 of pA_ligand
show sticks, pA_pocket
color lightorange, pA_pocket
set stick_radius, 0.15, pA_pocket
# Pocket surface
show surface, pA_pocket
color white, pA_pocket
set transparency, 0.40, pA_pocket

# Ligand
hide everything, pA_ligand
show sticks, pA_ligand
color red, pA_ligand
set stick_radius, 0.28, pA_ligand

# H-bonds
distance pA_hb, pA_receptor, pA_ligand, 3.5, mode=2
color yellow, pA_hb
hide labels, pA_hb

# Label
pseudoatom pA_label, pos=[0, 30, 0]
set label_size, 13
label pA_label, "(a) Hydromorphone  DG = -8.58"

# ═══ PANEL B: Biotin (−6.76) — Natural Substrate ═══
create pB_receptor, SMVT
create pB_ligand, Biotin
translate [80, 0, 0], pB_receptor
translate [80, 0, 0], pB_ligand

hide everything, pB_receptor
show cartoon, pB_receptor
color grey85, pB_receptor

select pB_pocket, pB_receptor within 5 of pB_ligand
show sticks, pB_pocket
color lightorange, pB_pocket
set stick_radius, 0.15, pB_pocket
show surface, pB_pocket
color white, pB_pocket
set transparency, 0.40, pB_pocket

hide everything, pB_ligand
show sticks, pB_ligand
color green, pB_ligand
set stick_radius, 0.28, pB_ligand

distance pB_hb, pB_receptor, pB_ligand, 3.5, mode=2
color yellow, pB_hb
hide labels, pB_hb

pseudoatom pB_label, pos=[80, 30, 0]
set label_size, 13
label pB_label, "(b) Biotin  DG = -6.76"

# ═══ PANEL C: All 8 Elite Hits — Chemical Diversity ═══
create pC_receptor, SMVT
translate [160, 0, 0], pC_receptor

hide everything, pC_receptor
show cartoon, pC_receptor
color grey85, pC_receptor
set cartoon_transparency, 0.4, pC_receptor

# Pocket
select pC_pocket, pC_receptor and resid 116-142+228-265+309-331+350-378
show surface, pC_pocket
color grey90, pC_pocket
set transparency, 0.55, pC_pocket

# All 8 elite hits in same pocket
python
import pymol
hits = ['HM','FUR','NAF','PHB','PTB','DIC','CAR','BTB']
colors = ['red','orange','marine','purple','violetpurple','yelloworange','deeppurple','magenta']
offset = [160, 0, 0]
for name, col in zip(hits, colors):
    pymol.cmd.create(f'pC_{name}', name)
    pymol.cmd.translate(offset, f'pC_{name}')
    pymol.cmd.hide('everything', f'pC_{name}')
    pymol.cmd.show('sticks', f'pC_{name}')
    pymol.cmd.color(col, f'pC_{name}')
    pymol.cmd.set('stick_radius', 0.22, f'pC_{name}')
python end

pseudoatom pC_label, pos=[160, 30, 0]
label pC_label, "(c) 8 Elite Hits Overlay"

# ═══ PANEL D: Controls Comparison ═══
create pD_receptor, SMVT
translate [240, 0, 0], pD_receptor

hide everything, pD_receptor
show cartoon, pD_receptor
color grey85, pD_receptor
set cartoon_transparency, 0.4, pD_receptor

select pD_pocket, pD_receptor and resid 116-142+228-265+309-331+350-378
show surface, pD_pocket
color grey90, pD_pocket
set transparency, 0.55, pD_pocket

# Controls
python
ctrls = ['Biotin','GBP','RIB','ESK']
colors = ['green','teal','grey70','cyan']
labels = ['Biotin (Ref)','Gabapentin (FDA+)','Riboflavin (-Ctrl)','Esketamine (Pilot)']
offset = [240, 0, 0]
for name in ctrls:
    pymol.cmd.create(f'pD_{name}', name)
    pymol.cmd.translate(offset, f'pD_{name}')
    pymol.cmd.hide('everything', f'pD_{name}')
    pymol.cmd.show('sticks', f'pD_{name}')
python end
color green, pD_Biotin
color teal, pD_GBP
color grey70, pD_RIB
color cyan, pD_ESK
set stick_radius, 0.22, pD_Biotin
set stick_radius, 0.22, pD_GBP
set stick_radius, 0.22, pD_RIB
set stick_radius, 0.22, pD_ESK

pseudoatom pD_label, pos=[240, 30, 0]
label pD_label, "(d) Controls + Pilot"

# ═══ Hide originals, zoom to all 4 panels ═══
hide everything, SMVT
hide everything, HM
hide everything, Biotin
hide everything, FUR
hide everything, NAF
hide everything, PHB
hide everything, PTB
hide everything, DIC
hide everything, CAR
hide everything, BTB
hide everything, GBP
hide everything, RIB
hide everything, ESK

zoom visible
extend visible, 1.1
move z, 40

# ═══ Legend ═══
print "╔══════════════════════════════════════════════════╗"
print "║  SMVT Docking — Publication Figure Layout      ║"
print "║  (a) Best Hit    (b) Natural Substrate          ║"
print "║  (c) 8 Elite Hits    (d) Controls               ║"
print "╠══════════════════════════════════════════════════╣"
print "║  ray 2400,1600; png fig_docking.png             ║"
print "╚══════════════════════════════════════════════════╝"
