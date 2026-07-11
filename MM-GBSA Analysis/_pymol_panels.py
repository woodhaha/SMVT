"""
PyMOL 4-panel binding pocket view → two separate figures (4 each).
"""
import sys, os, warnings
warnings.filterwarnings("ignore")

ANALYSIS = r'D:\Researching\SMVT\MM-GBSA Analysis'
ALL_COMPS = ['BIOTIN', 'HYDROMORPHONE', 'GABAPENTIN_ENACARBIL', 'NAFTAZONE',
             'ESKETAMINE', 'FUROSEMIDE', 'PHENOBARBITAL', 'RIBOFLAVIN']
LABELS = {'BIOTIN': 'Biotin (REF)', 'HYDROMORPHONE': 'Hydromorphone',
          'GABAPENTIN_ENACARBIL': 'Gabapentin Enacarbil', 'NAFTAZONE': 'Nafazone',
          'ESKETAMINE': 'Esketamine', 'FUROSEMIDE': 'Furosemide',
          'PHENOBARBITAL': 'Phenobarbital', 'RIBOFLAVIN': 'Riboflavin'}
COLORS = {'BIOTIN': 'blue', 'HYDROMORPHONE': 'magenta',
          'GABAPENTIN_ENACARBIL': 'green', 'NAFTAZONE': 'orange',
          'ESKETAMINE': 'red', 'FUROSEMIDE': 'yellow',
          'PHENOBARBITAL': 'purple', 'RIBOFLAVIN': 'salmon'}
SCORES = {'BIOTIN':(-6.38,-29.97), 'HYDROMORPHONE':(-7.72,-29.89),
          'GABAPENTIN_ENACARBIL':(-6.97,-43.33), 'NAFTAZONE':(-8.03,-22.36),
          'ESKETAMINE':(-7.36,-27.91), 'FUROSEMIDE':(-7.56,6.40),
          'PHENOBARBITAL':(-7.31,-23.61), 'RIBOFLAVIN':(-7.56,-41.48)}
TITTL = 'SMVT — Vina Docking Binding Modes'

# Split into 2 groups of 4
GROUPS = [ALL_COMPS[:4], ALL_COMPS[4:]]

import pymol
pymol.finish_launching(['pymol', '-cq'])
cmd = pymol.cmd
cmd.bg_color('white')
cmd.set('ray_opaque_background', 1)

cmd.load(os.path.join(ANALYSIS, 'BIOTIN_receptor.pdb'), 'receptor')
cmd.hide('everything', 'receptor')
cmd.show('cartoon', 'receptor')
cmd.color('gray80', 'receptor')
cmd.set('cartoon_transparency', 0.2, 'receptor')
cmd.select('pocket', 'receptor and chain A and resi 76+79+80+84+88+99+106+156+259+263+266+267+270+271+277+301+302+305+424+428+431+528')
cmd.show('sticks', 'pocket')
cmd.color('cyan', 'pocket')
cmd.set('stick_radius', 0.18)
cmd.orient('receptor')
cmd.zoom('receptor', -6)
common_view = cmd.get_view()

for name in ALL_COMPS:
    fn = os.path.join(ANALYSIS, f'{name}_docked.pdb')
    if not os.path.exists(fn):
        print(f'MISSING: {fn}, skipping')
        continue
    cmd.load(fn, 'lig')
    cmd.hide('everything', 'lig')
    cmd.show('sticks', 'lig')
    cmd.color(COLORS[name], 'lig')
    cmd.set('stick_radius', 0.4, 'lig')
    cmd.set_view(common_view)
    cmd.png(os.path.join(ANALYSIS, f'_panel_{name}.png'), width=3200, height=2800, dpi=300, ray=1)
    print(f'Rendered: {name}')
    cmd.delete('lig')
    cmd.delete('hb')

cmd.quit()

# ── PIL composite: two separate 2×2 figures ──
from PIL import Image, ImageDraw, ImageFont as PILFont

try:
    ft = PILFont.truetype('arial.ttf', 350)
    fl = PILFont.truetype('arial.ttf', 125)
except:
    ft = PILFont.truetype('arial.ttf', 100)
    fl = PILFont.truetype('arial.ttf', 50)

gap, title_h, bar_h = 80, 350, 200
pw, ph = 3200, 2800  # known panel size

for g_idx, comps in enumerate(GROUPS):
    W = pw * 2 + gap
    H = ph * 2 + gap + title_h
    img = Image.new('RGBA', (W, H), (255, 255, 255, 255))
    draw = ImageDraw.Draw(img)

    for i, name in enumerate(comps):
        panel = Image.open(os.path.join(ANALYSIS, f'_panel_{name}.png'))
        col = i % 2; row = i // 2
        x = col * pw + (gap // 2 if col > 0 else 0)
        y = row * ph + gap * row + title_h
        img.paste(panel, (x, y))

        bx = x; by = y + ph - bar_h
        v, m = SCORES[name]
        text = f'{LABELS[name]}    Vina={v:.2f}    MMGBSA={m:.1f}'
        overlay = Image.new('RGBA', (pw, bar_h), (255, 255, 255, 220))
        img.paste(overlay, (bx, by), overlay)
        draw.text((bx + 40, by + 30), text, fill='#000', font=fl)

    draw.text((40, 15), TITTL, fill='black', font=ft)

    out_fn = f'ALL_pymol_4panel_g{g_idx+1}.png'
    out = os.path.join(ANALYSIS, out_fn)
    img.save(out, dpi=(400, 400))
    print(f'Saved: {out}')

for c in ALL_COMPS:
    os.remove(os.path.join(ANALYSIS, f'_panel_{c}.png'))
print('Temp cleaned')
