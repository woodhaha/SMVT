cd D:/Researching/SMVT/03_Analysis/pymol
bg_color white
set valence, off
set grid_mode, on

# ── Load receptor (reused) ──
load SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey80, SMVT

# ── Load ligands, translate apart ──
load Biotin_docked.pdb, Biotin
load Phenobarbital_docked.pdb, Phenobarbital
load Naftazone_docked.pdb, Naftazone
load Furosemide_docked.pdb, Furosemide
load Diclofenac_docked.pdb, Diclofenac
load Esketamine_docked.pdb, Esketamine

# Translate each to its own position (grid layout)
translate [0, 0, 0], Biotin
translate [0, 25, 0], Phenobarbital
translate [0, 50, 0], Naftazone
translate [0, 75, 0], Furosemide
translate [0, 100, 0], Diclofenac
translate [0, 125, 0], Esketamine

# Copy SMVT to each ligand position
create SMVT_biotin, SMVT
create SMVT_pheno, SMVT
create SMVT_naft, SMVT
create SMVT_furo, SMVT
create SMVT_dicl, SMVT
create SMVT_esket, SMVT

translate [0, 0, 0], SMVT_biotin
translate [0, 25, 0], SMVT_pheno
translate [0, 50, 0], SMVT_naft
translate [0, 75, 0], SMVT_furo
translate [0, 100, 0], SMVT_dicl
translate [0, 125, 0], SMVT_esket

# Style ligands
hide everything, Biotin
show sticks, Biotin
color green, Biotin
set stick_radius, 0.25, Biotin

hide everything, Phenobarbital
show sticks, Phenobarbital
color purple, Phenobarbital
set stick_radius, 0.25, Phenobarbital

hide everything, Naftazone
show sticks, Naftazone
color blue, Naftazone
set stick_radius, 0.25, Naftazone

hide everything, Furosemide
show sticks, Furosemide
color red, Furosemide
set stick_radius, 0.25, Furosemide

hide everything, Diclofenac
show sticks, Diclofenac
color orange, Diclofenac
set stick_radius, 0.25, Diclofenac

hide everything, Esketamine
show sticks, Esketamine
color cyan, Esketamine
set stick_radius, 0.25, Esketamine

# Labels
label SMVT_biotin and name CA and resi 200, "BIOTIN (-6.76)"
set label_color, green
label SMVT_pheno and name CA and resi 200, "PHENOBARBITAL (-8.30)"
set label_color, purple
label SMVT_naft and name CA and resi 200, "NAFTAZONE (-8.34)"
set label_color, blue
label SMVT_furo and name CA and resi 200, "FUROSEMIDE (-8.36)"
set label_color, red
label SMVT_dicl and name CA and resi 200, "DICLOFENAC (-7.15)"
set label_color, orange
label SMVT_esket and name CA and resi 200, "ESKETAMINE (-7.58)"
set label_color, cyan

# Hide original SMVT
hide everything, SMVT

# Zoom to see all
zoom visible

# Print legend
print "============================================"
print "  SMVT Docking — Split View"
print "============================================"
print "  GREEN:   Biotin        (-6.76)  Natural substrate"
print "  PURPLE:  Phenobarbital (-8.30)  Barbiturate 100% hit"
print "  BLUE:    Naftazone     (-8.34)  Naphthoquinone"
print "  RED:     Furosemide    (-8.36)  Sulfonamide"
print "  ORANGE:  Diclofenac    (-7.15)  NSAID"
print "  CYAN:    Esketamine    (-7.58)  Arylcyclohexylamine"
print "============================================"
