# ═══════════════════════════════════════════════════════
# SMVT (SLC5A6) Docking — All 8 Elite Hits + Controls
# ═══════════════════════════════════════════════════════
# Run: pymol D:\Researching\SMVT\03_Analysis\pymol\smvt_all_hits.pml

cd D:/Researching/SMVT/03_Analysis/pymol
bg_color white
set valence, off
set grid_mode, on

# ═══ Receptor ═══
load SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey80, SMVT

# ── Highlight binding pocket residues (transmembrane domain) ──
select pocket, SMVT and resid 116-142+228-265+309-331+350-378
show surface, pocket
color grey90, pocket
set transparency, 0.65, pocket
select pocket_cartoon, SMVT and resid 116-142+228-265+309-331+350-378
color salmon, pocket_cartoon

# ═══ Load all 12 ligands ═══
load Hydromorphone_docked.pdb, Hydromorphone
load Furosemide_docked.pdb, Furosemide
load Naftazone_docked.pdb, Naftazone
load Phenobarbital_docked.pdb, Phenobarbital
load Pentobarbital_docked.pdb, Pentobarbital
load Diclofenac_docked.pdb, Diclofenac
load Carprofen_docked.pdb, Carprofen
load Butabarbital_docked.pdb, Butabarbital
load Biotin_docked.pdb, Biotin
load Gabapentin_enacarbil_docked.pdb, Gabapentin
load Riboflavin_docked.pdb, Riboflavin
load Esketamine_docked.pdb, Esketamine

# ═══ Style all ligands ═══
# Helper macro
macro style_ligand, ligand, col, rad {
    hide everything, ligand
    show sticks, ligand
    color col, ligand
    set stick_radius, rad, ligand
}

# Top 8 Hits — warm colors (best affinity → weaker)
style_ligand Hydromorphone, deeppink, 0.25      # −8.58 #1
style_ligand Furosemide, red, 0.25               # −8.36 #2
style_ligand Naftazone, marine, 0.25             # −8.34 #3
style_ligand Phenobarbital, purple, 0.25         # −8.30 #4
style_ligand Pentobarbital, violetpurple, 0.25   # −8.18 #5
style_ligand Diclofenac, orange, 0.25            # −8.07 #6
style_ligand Carprofen, yelloworange, 0.25       # −8.04 #7
style_ligand Butabarbital, magenta, 0.25         # −8.02 #8

# Controls — cool colors
style_ligand Biotin, green, 0.25                 # Natural substrate
style_ligand Gabapentin, teal, 0.25              # FDA positive control
style_ligand Riboflavin, grey60, 0.25            # Negative control
style_ligand Esketamine, cyan, 0.25              # Pilot compound

# ═══ Groups ═══
group Elite_Hits, Hydromorphone Furosemide Naftazone Phenobarbital Pentobarbital Diclofenac Carprofen Butabarbital
group Controls, Biotin Gabapentin Riboflavin
group Pilot, Esketamine

# ═══ Label Binding Site ═══
label SMVT and name CA and resi 150, "SMVT Binding Pocket"
set label_color, grey40
set label_size, 24

# ═══ H-bond detection (show potential interactions) ═══
# Biotin ureido in pocket
distance hb_biotin, SMVT, Biotin, 3.5, mode=2
color green, hb_biotin
hide labels, hb_biotin

# ═══ View setup ═══
zoom SMVT
extend SMVT, 4.0
set antialias, 2
set depth_cue, on
set ray_trace_mode, 1

# ═══ Print legend ═══
print "╔══════════════════════════════════════════════╗"
print "║  SMVT Docking — 8 Elite Hits + 1 Control    ║"
print "╠══════════════════════════════════════════════╣"
print "║ ELITE HITS (dG < -8.0 kcal/mol):            ║"
print "║  PINK:   Hydromorphone   -8.58 #1 opioid    ║"
print "║  RED:    Furosemide      -8.36 sulfonamide  ║"
print "║  BLUE:   Naftazone       -8.34 naphthoquin  ║"
print "║  PURPLE: Phenobarbital   -8.30 barbiturate  ║"
print "║  VIOLET: Pentobarbital   -8.18 barbiturate  ║"
print "║  ORANGE: Diclofenac      -8.07 NSAID        ║"
print "║  GOLD:   Carprofen       -8.04 NSAID        ║"
print "║  MAGENTA:Butabarbital    -8.02 barbiturate  ║"
print "╠══════════════════════════════════════════════╣"
print "║ CONTROLS:                                   ║"
print "║  GREEN:  Biotin (-6.76) natural substrate   ║"
print "║  TEAL:   Gabapentin enacarbil (-6.63)       ║"
print "║  GREY:   Riboflavin (-0.01) non-binder      ║"
print "╚══════════════════════════════════════════════╝"
print ""
print "Commands:"
print "  disable Elite_Hits       # hide all hits"
print "  enable Elite_Hits        # show all hits"
print "  disable Controls         # hide controls"
print "  zoom Hydromorphone       # focus on best hit"
print "  ray 2400,1600; png docking_all.png   # render"
print "  save smvt_docking.pse    # save session"
