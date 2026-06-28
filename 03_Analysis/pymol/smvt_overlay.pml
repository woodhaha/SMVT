# ═══════════════════════════════════════════════════════
# SMVT — 8 Ligands Overlay in Single Binding Pocket
# ═══════════════════════════════════════════════════════

reinitialize
bg_color white
set valence, off
set antialias, 2
set depth_cue, on
set cartoon_fancy_helices, on
set ray_trace_mode, 1

# ═══ Receptor ═══
load D:/Researching/SMVT/03_Analysis/pymol/SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey85, SMVT
set cartoon_transparency, 0.3, SMVT

# ── Pocket surface ──
select pocket, SMVT and resid 116-142+228-265+309-331+350-378
show surface, pocket
color grey90, pocket
set transparency, 0.50, pocket

# ═══ 8 Ligands ═══
load D:/Researching/SMVT/03_Analysis/pymol/Hydromorphone_docked.pdb, Hydromorphone
load D:/Researching/SMVT/03_Analysis/pymol/Furosemide_docked.pdb, Furosemide
load D:/Researching/SMVT/03_Analysis/pymol/Naftazone_docked.pdb, Naftazone
load D:/Researching/SMVT/03_Analysis/pymol/Phenobarbital_docked.pdb, Phenobarbital
load D:/Researching/SMVT/03_Analysis/pymol/Biotin_docked.pdb, Biotin
load D:/Researching/SMVT/03_Analysis/pymol/Gabapentin_enacarbil_docked.pdb, Gabapentin
load D:/Researching/SMVT/03_Analysis/pymol/Riboflavin_docked.pdb, Riboflavin
load D:/Researching/SMVT/03_Analysis/pymol/Esketamine_docked.pdb, Esketamine

# ═══ Style — Test Compounds ═══
hide everything, Hydromorphone
show sticks, Hydromorphone
color red, Hydromorphone
set stick_radius, 0.30, Hydromorphone

hide everything, Furosemide
show sticks, Furosemide
color orange, Furosemide
set stick_radius, 0.30, Furosemide

hide everything, Naftazone
show sticks, Naftazone
color marine, Naftazone
set stick_radius, 0.30, Naftazone

hide everything, Phenobarbital
show sticks, Phenobarbital
color purple, Phenobarbital
set stick_radius, 0.30, Phenobarbital

# ═══ Style — Controls ═══
hide everything, Biotin
show sticks, Biotin
color green, Biotin
set stick_radius, 0.28, Biotin

hide everything, Gabapentin
show sticks, Gabapentin
color teal, Gabapentin
set stick_radius, 0.28, Gabapentin

hide everything, Riboflavin
show sticks, Riboflavin
color grey70, Riboflavin
set stick_radius, 0.28, Riboflavin

# ═══ Style — Pilot ═══
hide everything, Esketamine
show sticks, Esketamine
color cyan, Esketamine
set stick_radius, 0.28, Esketamine

# ═══ H-bond detection (Biotin as reference) ═══
distance hb_biotin, SMVT, Biotin, 3.5, mode=2
color yellow, hb_biotin
hide labels, hb_biotin

# ═══ Groups ═══
group TestCompounds, Hydromorphone Furosemide Naftazone Phenobarbital
group Controls, Biotin Gabapentin Riboflavin
group Pilot, Esketamine

# ═══ Zoom to pocket ═══
zoom pocket
extend pocket, 3.0

# ═══ Label ═══
set label_size, 20
pseudoatom title, pos=[0, 0, 0]
# hide the pseudoatom itself
hide nonbonded, title

# ═══ Scene preset ═══
# Front view of pocket
turn x, 20
turn y, -10

# ═══ Legend in console ═══
print "╔════════════════════════════════════════╗"
print "║  SMVT Binding Pocket — 8 Ligands     ║"
print "╠════════════════════════════════════════╣"
print "║ TEST COMPOUNDS (dG < -8.0):           ║"
print "║  RED:    Hydromorphone  -8.58         ║"
print "║  ORANGE: Furosemide     -8.36         ║"
print "║  BLUE:   Naftazone      -8.34         ║"
print "║  PURPLE: Phenobarbital  -8.30         ║"
print "╠════════════════════════════════════════╣"
print "║ CONTROLS:                             ║"
print "║  GREEN:  Biotin         -6.76 (Ref)   ║"
print "║  TEAL:   Gabapentin     -6.63 (FDA)   ║"
print "║  GREY:   Riboflavin     -0.01 (-Ctrl) ║"
print "╠════════════════════════════════════════╣"
print "║ PILOT:                                ║"
print "║  CYAN:   Esketamine     -7.58         ║"
print "╚════════════════════════════════════════╝"
print ""
print "Toggle groups:"
print "  disable TestCompounds   # hide test hits"
print "  disable Controls        # hide controls"
print "  enable TestCompounds    # show again"
print ""
print "Focus:"
print "  zoom Hydromorphone      # best hit"
print "  zoom Biotin             # natural substrate"
print ""
print "Export: ray 2400,1600; png smvt_overlay.png"
