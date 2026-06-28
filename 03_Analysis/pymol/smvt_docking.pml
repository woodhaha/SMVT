cd D:/Researching/SMVT/03_Analysis/pymol
bg_color white
set valence, off

load SMVT_receptor.pdb, SMVT
hide everything, SMVT
show cartoon, SMVT
color grey80, SMVT

load Biotin_docked.pdb, Biotin
load Phenobarbital_docked.pdb, Phenobarbital
load Naftazone_docked.pdb, Naftazone
load Furosemide_docked.pdb, Furosemide
load Diclofenac_docked.pdb, Diclofenac
load Esketamine_docked.pdb, Esketamine

hide everything, Biotin
show sticks, Biotin
color green, Biotin
set stick_radius, 0.2, Biotin

hide everything, Phenobarbital
show sticks, Phenobarbital
color purple, Phenobarbital
set stick_radius, 0.2, Phenobarbital

hide everything, Naftazone
show sticks, Naftazone
color blue, Naftazone
set stick_radius, 0.2, Naftazone

hide everything, Furosemide
show sticks, Furosemide
color red, Furosemide
set stick_radius, 0.2, Furosemide

hide everything, Diclofenac
show sticks, Diclofenac
color orange, Diclofenac
set stick_radius, 0.2, Diclofenac

hide everything, Esketamine
show sticks, Esketamine
color cyan, Esketamine
set stick_radius, 0.2, Esketamine

group Hits, Furosemide Naftazone Phenobarbital Diclofenac Esketamine
group Controls, Biotin

zoom SMVT
orient SMVT
