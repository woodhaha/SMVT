# ═══════════════════════════════════════════════════════
# SMVT Docking — Export Publication Panels
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
set label_size, 12

# ═══ Load assets ═══
load D:/Researching/SMVT/03_Analysis/pymol/SMVT_receptor.pdb, SMVT
load D:/Researching/SMVT/03_Analysis/pymol/Hydromorphone_docked.pdb, HM
load D:/Researching/SMVT/03_Analysis/pymol/Biotin_docked.pdb, Biotin
load D:/Researching/SMVT/03_Analysis/pymol/Furosemide_docked.pdb, FUR
load D:/Researching/SMVT/03_Analysis/pymol/Naftazone_docked.pdb, NAF
load D:/Researching/SMVT/03_Analysis/pymol/Phenobarbital_docked.pdb, PHB
load D:/Researching/SMVT/03_Analysis/pymol/Pentobarbital_docked.pdb, PTB
load D:/Researching/SMVT/03_Analysis/pymol/Diclofenac_docked.pdb, DIC
load D:/Researching/SMVT/03_Analysis/pymol/Carprofen_docked.pdb, CAR
load D:/Researching/SMVT/03_Analysis/pymol/Butabarbital_docked.pdb, BTB
load D:/Researching/SMVT/03_Analysis/pymol/Gabapentin_enacarbil_docked.pdb, GBP
load D:/Researching/SMVT/03_Analysis/pymol/Riboflavin_docked.pdb, RIB
load D:/Researching/SMVT/03_Analysis/pymol/Esketamine_docked.pdb, ESK

# ═══════════════════════════════════════
# PANEL A: Hydromorphone (best hit)
# ═══════════════════════════════════════
create pA_r, SMVT
create pA_l, HM

hide everything, pA_r
show cartoon, pA_r
color grey85, pA_r

select pA_pocket, pA_r within 5 of pA_l
show sticks, pA_pocket
color lightorange, pA_pocket
set stick_radius, 0.14, pA_pocket
show surface, pA_pocket
color white, pA_pocket
set transparency, 0.40, pA_pocket

hide everything, pA_l
show sticks, pA_l
color red, pA_l
set stick_radius, 0.28, pA_l

distance pA_hb, pA_r, pA_l, 3.5, mode=2
color yellow, pA_hb
hide labels, pA_hb

pseudoatom pA_lbl, pos=[0, 28, 0]
label pA_lbl, "(a) Hydromorphone  DG = -8.58"

# ═══════════════════════════════════════
# PANEL B: Biotin (natural substrate)
# ═══════════════════════════════════════
create pB_r, SMVT
create pB_l, Biotin
translate [80, 0, 0], pB_r
translate [80, 0, 0], pB_l

hide everything, pB_r
show cartoon, pB_r
color grey85, pB_r

select pB_pocket, pB_r within 5 of pB_l
show sticks, pB_pocket
color lightorange, pB_pocket
set stick_radius, 0.14, pB_pocket
show surface, pB_pocket
color white, pB_pocket
set transparency, 0.40, pB_pocket

hide everything, pB_l
show sticks, pB_l
color green, pB_l
set stick_radius, 0.28, pB_l

distance pB_hb, pB_r, pB_l, 3.5, mode=2
color yellow, pB_hb
hide labels, pB_hb

pseudoatom pB_lbl, pos=[80, 28, 0]
label pB_lbl, "(b) Biotin  DG = -6.76"

# ═══════════════════════════════════════
# PANEL C: 8 Elite Hits Overlay
# ═══════════════════════════════════════
create pC_r, SMVT
translate [160, 0, 0], pC_r

hide everything, pC_r
show cartoon, pC_r
color grey85, pC_r
set cartoon_transparency, 0.35, pC_r

select pC_pocket, pC_r and resid 116-142+228-265+309-331+350-378
show surface, pC_pocket
color grey90, pC_pocket
set transparency, 0.55, pC_pocket

# Clone 8 hits to panel C
create pC_HM, HM
create pC_FUR, FUR
create pC_NAF, NAF
create pC_PHB, PHB
create pC_PTB, PTB
create pC_DIC, DIC
create pC_CAR, CAR
create pC_BTB, BTB

translate [160, 0, 0], pC_HM
translate [160, 0, 0], pC_FUR
translate [160, 0, 0], pC_NAF
translate [160, 0, 0], pC_PHB
translate [160, 0, 0], pC_PTB
translate [160, 0, 0], pC_DIC
translate [160, 0, 0], pC_CAR
translate [160, 0, 0], pC_BTB

hide everything, pC_HM
show sticks, pC_HM
color red, pC_HM
set stick_radius, 0.22, pC_HM

hide everything, pC_FUR
show sticks, pC_FUR
color orange, pC_FUR
set stick_radius, 0.22, pC_FUR

hide everything, pC_NAF
show sticks, pC_NAF
color blue, pC_NAF
set stick_radius, 0.22, pC_NAF

hide everything, pC_PHB
show sticks, pC_PHB
color purple, pC_PHB
set stick_radius, 0.22, pC_PHB

hide everything, pC_PTB
show sticks, pC_PTB
color violetpurple, pC_PTB
set stick_radius, 0.22, pC_PTB

hide everything, pC_DIC
show sticks, pC_DIC
color yelloworange, pC_DIC
set stick_radius, 0.22, pC_DIC

hide everything, pC_CAR
show sticks, pC_CAR
color deeppurple, pC_CAR
set stick_radius, 0.22, pC_CAR

hide everything, pC_BTB
show sticks, pC_BTB
color magenta, pC_BTB
set stick_radius, 0.22, pC_BTB

pseudoatom pC_lbl, pos=[160, 28, 0]
label pC_lbl, "(c) 8 Elite Hits Overlay  (dG < -8.0)"

# ═══════════════════════════════════════
# PANEL D: Controls + Pilot
# ═══════════════════════════════════════
create pD_r, SMVT
translate [240, 0, 0], pD_r

hide everything, pD_r
show cartoon, pD_r
color grey85, pD_r
set cartoon_transparency, 0.35, pD_r

select pD_pocket, pD_r and resid 116-142+228-265+309-331+350-378
show surface, pD_pocket
color grey90, pD_pocket
set transparency, 0.55, pD_pocket

# Clone controls to panel D
create pD_Bio, Biotin
create pD_GBP, GBP
create pD_RIB, RIB
create pD_ESK, ESK

translate [240, 0, 0], pD_Bio
translate [240, 0, 0], pD_GBP
translate [240, 0, 0], pD_RIB
translate [240, 0, 0], pD_ESK

hide everything, pD_Bio
show sticks, pD_Bio
color green, pD_Bio
set stick_radius, 0.22, pD_Bio

hide everything, pD_GBP
show sticks, pD_GBP
color teal, pD_GBP
set stick_radius, 0.22, pD_GBP

hide everything, pD_RIB
show sticks, pD_RIB
color grey70, pD_RIB
set stick_radius, 0.22, pD_RIB

hide everything, pD_ESK
show sticks, pD_ESK
color cyan, pD_ESK
set stick_radius, 0.22, pD_ESK

pseudoatom pD_lbl, pos=[240, 28, 0]
label pD_lbl, "(d) Controls + Pilot"

# ═══ Hide originals ═══
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

# ═══ Zoom ═══
zoom visible
extend visible, 1.1
move z, 45

# ═══ RENDER ═══
print "Rendering panels..."

# Full 4-panel
viewport 2400, 1200
ray 2400, 1200
png D:/Researching/SMVT/03_Analysis/pymol/exports/Fig_docking_4panel.png
print "[1/5] Fig_docking_4panel.png"

# Panel A
zoom pA_r
extend pA_r, 2.0
viewport 1200, 1200
ray 1200, 1200
png D:/Researching/SMVT/03_Analysis/pymol/exports/Panel_A_Hydromorphone.png
print "[2/5] Panel_A_Hydromorphone.png"

# Panel B
zoom pB_r
extend pB_r, 2.0
ray 1200, 1200
png D:/Researching/SMVT/03_Analysis/pymol/exports/Panel_B_Biotin.png
print "[3/5] Panel_B_Biotin.png"

# Panel C
zoom pC_r
extend pC_r, 2.0
ray 1200, 1200
png D:/Researching/SMVT/03_Analysis/pymol/exports/Panel_C_8Hits.png
print "[4/5] Panel_C_8Hits.png"

# Panel D
zoom pD_r
extend pD_r, 2.0
ray 1200, 1200
png D:/Researching/SMVT/03_Analysis/pymol/exports/Panel_D_Controls.png
print "[5/5] Panel_D_Controls.png"

print "DONE: D:/Researching/SMVT/03_Analysis/pymol/exports/"
