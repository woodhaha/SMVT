"""
SMVT Docking — Render publication panels via PyMOL Python API.
Run: python render_export.py
"""
import os, sys
os.chdir(os.path.join(os.path.dirname(__file__), "..", ".."))
os.makedirs("exports", exist_ok=True)

# Launch PyMOL in headless mode (no GUI, faster rendering)
import pymol
pymol.pymol_argv = ['pymol', '-qc']  # -q = quiet, -c = command line (no GUI)
pymol.finish_launching()

cmd = pymol.cmd

# ── Settings ──
cmd.bg_color('white')
cmd.set('valence', 'off')
cmd.set('antialias', 3)
cmd.set('ray_trace_mode', 1)
cmd.set('ray_trace_gain', 0)
cmd.set('cartoon_fancy_helices', 'on')
cmd.set('stick_quality', 15)
cmd.set('label_size', 12)

# ── Load ──
cmd.load('SMVT_receptor.pdb', 'SMVT')
cmd.load('Hydromorphone_docked.pdb', 'HM')
cmd.load('Biotin_docked.pdb', 'Biotin')
cmd.load('Furosemide_docked.pdb', 'FUR')
cmd.load('Naftazone_docked.pdb', 'NAF')
cmd.load('Phenobarbital_docked.pdb', 'PHB')
cmd.load('Pentobarbital_docked.pdb', 'PTB')
cmd.load('Diclofenac_docked.pdb', 'DIC')
cmd.load('Carprofen_docked.pdb', 'CAR')
cmd.load('Butabarbital_docked.pdb', 'BTB')
cmd.load('Gabapentin_enacarbil_docked.pdb', 'GBP')
cmd.load('Riboflavin_docked.pdb', 'RIB')
cmd.load('Esketamine_docked.pdb', 'ESK')
print(f"Loaded {len(cmd.get_object_list())} objects: {cmd.get_object_list()}")

def style_receptor(rec, lig):
    """Style receptor with pocket highlight around ligand"""
    cmd.hide('everything', rec)
    cmd.show('cartoon', rec)
    cmd.color('grey85', rec)
    sel = f'{rec} within 5 of {lig}'
    cmd.select('tmp_pocket', sel)
    cmd.show('sticks', 'tmp_pocket')
    cmd.color('lightorange', 'tmp_pocket')
    cmd.set('stick_radius', 0.14, 'tmp_pocket')
    cmd.show('surface', 'tmp_pocket')
    cmd.color('white', 'tmp_pocket')
    cmd.set('transparency', 0.40, 'tmp_pocket')
    # H-bonds
    cmd.distance('tmp_hb', rec, lig, 3.5, mode=2)
    cmd.color('yellow', 'tmp_hb')
    cmd.hide('labels', 'tmp_hb')

def style_ligand(lig, color):
    cmd.hide('everything', lig)
    cmd.show('sticks', lig)
    cmd.color(color, lig)
    cmd.set('stick_radius', 0.28, lig)

def add_label(rec, lig, text):
    cmd.pseudoatom(f'{rec}_lbl', pos=[0, 28, 0])
    cmd.label(f'{rec}_lbl', f'"{text}"')

# ═══ PANEL A ═══
print("Panel A: Hydromorphone")
cmd.create('pA_r', 'SMVT')
cmd.create('pA_l', 'HM')
style_receptor('pA_r', 'pA_l')
style_ligand('pA_l', 'red')
add_label('pA_r', 'pA_l', "(a) Hydromorphone  DG = -8.58")

# ═══ PANEL B ═══
print("Panel B: Biotin")
cmd.create('pB_r', 'SMVT')
cmd.create('pB_l', 'Biotin')
cmd.translate([80, 0, 0], 'pB_r')
cmd.translate([80, 0, 0], 'pB_l')
style_receptor('pB_r', 'pB_l')
style_ligand('pB_l', 'green')
add_label('pB_r', 'pB_l', "(b) Biotin  DG = -6.76")

# ═══ PANEL C: 8 Hits ═══
print("Panel C: 8 Elite Hits")
cmd.create('pC_r', 'SMVT')
cmd.translate([160, 0, 0], 'pC_r')
cmd.hide('everything', 'pC_r')
cmd.show('cartoon', 'pC_r')
cmd.color('grey85', 'pC_r')
cmd.set('cartoon_transparency', 0.35, 'pC_r')
cmd.select('tmp_pkt', 'pC_r and resid 116-142+228-265+309-331+350-378')
cmd.show('surface', 'tmp_pkt')
cmd.color('grey90', 'tmp_pkt')
cmd.set('transparency', 0.55, 'tmp_pkt')

hits = [('HM','red'),('FUR','orange'),('NAF','blue'),('PHB','purple'),
        ('PTB','violetpurple'),('DIC','yelloworange'),('CAR','deeppurple'),('BTB','magenta')]
for name, col in hits:
    obj = f'pC_{name}'
    cmd.create(obj, name)
    cmd.translate([160, 0, 0], obj)
    cmd.hide('everything', obj)
    cmd.show('sticks', obj)
    cmd.color(col, obj)
    cmd.set('stick_radius', 0.22, obj)

cmd.pseudoatom('pC_lbl', pos=[160, 28, 0])
cmd.label('pC_lbl', '"(c) 8 Elite Hits  (dG < -8.0)"')

# ═══ PANEL D: Controls ═══
print("Panel D: Controls")
cmd.create('pD_r', 'SMVT')
cmd.translate([240, 0, 0], 'pD_r')
cmd.hide('everything', 'pD_r')
cmd.show('cartoon', 'pD_r')
cmd.color('grey85', 'pD_r')
cmd.set('cartoon_transparency', 0.35, 'pD_r')
cmd.select('tmp_pkt2', 'pD_r and resid 116-142+228-265+309-331+350-378')
cmd.show('surface', 'tmp_pkt2')
cmd.color('grey90', 'tmp_pkt2')
cmd.set('transparency', 0.55, 'tmp_pkt2')

ctrls = [('Biotin','green'),('GBP','teal'),('RIB','grey70'),('ESK','cyan')]
for name, col in ctrls:
    obj = f'pD_{name}'
    cmd.create(obj, name)
    cmd.translate([240, 0, 0], obj)
    cmd.hide('everything', obj)
    cmd.show('sticks', obj)
    cmd.color(col, obj)
    cmd.set('stick_radius', 0.22, obj)

cmd.pseudoatom('pD_lbl', pos=[240, 28, 0])
cmd.label('pD_lbl', '"(d) Controls + Pilot"')

# ═══ Hide originals ═══
for obj in cmd.get_object_list('all'):
    if not any(obj.startswith(p) for p in ['pA_','pB_','pC_','pD_']):
        cmd.hide('everything', obj)

# ═══ Render ═══
cmd.zoom('visible')
cmd.extend('visible', 1.1)
cmd.move('z', 45)

# Full 4-panel
print("Rendering 4-panel figure...")
cmd.viewport(2400, 1200)
cmd.ray(2400, 1200)
cmd.png('exports/Fig_docking_4panel.png')
print("  [1/5] Fig_docking_4panel.png ✓")

# Individual panels
panels = [
    ('pA_r', 'exports/Panel_A_Hydromorphone.png'),
    ('pB_r', 'exports/Panel_B_Biotin.png'),
    ('pC_r', 'exports/Panel_C_8Hits.png'),
    ('pD_r', 'exports/Panel_D_Controls.png'),
]
for i, (rec, out) in enumerate(panels, 2):
    print(f"  [{i}/5] {out}...")
    cmd.zoom(rec)
    cmd.extend(rec, 2.0)
    cmd.ray(1200, 1200)
    cmd.png(out)

print("\nDONE — all images in exports/")
print("  " + "\n  ".join(os.listdir('exports')))

cmd.quit()
