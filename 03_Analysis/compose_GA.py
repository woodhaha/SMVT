#!/usr/bin/env python3
"""SMVT 图形摘要合成 — 1328×531 (5:2 横版)。

四区布局（左→右）：
  1) 泛癌表达热图（读 pancancer_heatmap.png）
  2) 膜 + SMVT 转运体示意 + 生物素
  3) 虚拟筛选流程（3 步）
  4) 巴比妥酸 vs 生物素 ureido 环对比 + 机制轴

分子用 RDKit 2D 坐标 + matplotlib 手绘（无需 Cairo）。
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Ellipse, Circle
from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem

HEATMAP = r"D:\Researching\SMVT\03_Analysis\outputs\pancancer_heatmap.png"
OUT = r"D:\Researching\SMVT\03_Analysis\outputs\SMVT_graphical_abstract"

BIOTIN = "C1[C@H]2[C@@H]([C@@H](S1)CCCCC(=O)O)NC(=O)N2"
BARB = "O=C1NC(=O)CC(=O)N1"

ATOM_COLORS = {"C": "#37474f", "O": "#e53935", "N": "#1e88e5",
               "S": "#f9a825", "H": "#b0bec5"}

plt.rcParams.update({"font.family": "Arial", "font.size": 9})


def draw_molecule(ax, smiles, atom_r=0.16):
    """在 ax 里手绘分子骨架结构（RDKit 2D 坐标）。"""
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    conf = mol.GetConformer()
    pts = {i: (conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y)
           for i in range(mol.GetNumAtoms())}

    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        (x1, y1), (x2, y2) = pts[i], pts[j]
        order = bond.GetBondTypeAsDouble()
        if order > 1.0:  # 双键/芳香键：两条平行偏移线
            dx, dy = x2 - x1, y2 - y1
            L = (dx * dx + dy * dy) ** 0.5 or 1.0
            ox, oy = -dy / L * 0.06, dx / L * 0.06
            ax.plot([x1 + ox, x2 + ox], [y1 + oy, y2 + oy],
                    color="#37474f", lw=1.8, solid_capstyle="round", zorder=2)
            ax.plot([x1 - ox, x2 - ox], [y1 - oy, y2 - oy],
                    color="#37474f", lw=1.8, solid_capstyle="round", zorder=2)
        else:
            ax.plot([x1, x2], [y1, y2],
                    color="#37474f", lw=1.8, solid_capstyle="round", zorder=2)

    for i, atom in enumerate(mol.GetAtoms()):
        x, y = pts[i]
        sym = atom.GetSymbol()
        ax.add_patch(Circle((x, y), atom_r, color=ATOM_COLORS.get(sym, "#37474f"),
                            ec="white", lw=0.5, zorder=3))
        if sym != "C":  # 杂原子标注（碳省略）
            ax.text(x, y, sym, ha="center", va="center", fontsize=7.5,
                    color="white", fontweight="bold", zorder=4)

    xs = [p[0] for p in pts.values()]; ys = [p[1] for p in pts.values()]
    ax.set_xlim(min(xs) - 0.5, max(xs) + 0.5)
    ax.set_ylim(min(ys) - 0.5, max(ys) + 0.5)
    ax.set_aspect("equal")
    ax.axis("off")


# ── 画布 1328×531 ──
fig = plt.figure(figsize=(13.28, 5.31), dpi=100)

# 1) 泛癌热图（左，竖版）
ax_heat = fig.add_axes([0.010, 0.08, 0.155, 0.86])
ax_heat.imshow(plt.imread(HEATMAP))
ax_heat.axis("off")

# 2) 膜 + 转运体（左中）
ax_mem = fig.add_axes([0.185, 0.06, 0.235, 0.88])
ax_mem.set_xlim(0, 1); ax_mem.set_ylim(0, 1); ax_mem.axis("off")

for y in (0.46, 0.54):  # 脂双层
    ax_mem.plot([0.04, 0.96], [y, y], color="#9aa0a6", lw=2.5, zorder=2)
    for x in np.linspace(0.05, 0.95, 12):
        ax_mem.add_patch(Circle((x, y + (0.028 if y == 0.46 else -0.028)),
                                0.022, color="#cfd8dc", ec="#78909c", lw=0.8, zorder=3))

ax_mem.add_patch(Ellipse((0.5, 0.50), 0.62, 0.42, color="#e3f2fd",
                         ec="#0288d1", lw=2.5, zorder=4))

for i in range(13):  # 13 条 TM 螺旋
    ang = 2 * np.pi * i / 13
    ax_mem.add_patch(Circle((0.5 + 0.27 * np.cos(ang), 0.50 + 0.18 * np.sin(ang)),
                            0.035, color="#4fc3f7", ec="#0277bd", lw=1.2, zorder=5))

ax_mem.plot([0.44, 0.50, 0.56], [0.60, 0.68, 0.60], color="#01579b",
            lw=3, solid_capstyle="round", zorder=6)  # 扭结 TM1

for (nx, ny) in ((0.66, 0.62), (0.72, 0.55)):  # 2 Na+
    ax_mem.add_patch(Circle((nx, ny), 0.030, color="#7e57c2", ec="#4527a0",
                            lw=1, zorder=6))

ax_mem.text(0.5, 0.985, "SMVT / SLC5A6", ha="center", va="top",
            fontsize=11, fontweight="bold", color="#01579b")
ax_mem.text(0.50, 0.02, "biotin  ·  2 Na$^+$", ha="center", va="bottom",
            fontsize=8.5, color="#455a64")
ax_mem.text(0.90, 0.985, "kinked TM1", ha="right", va="top",
            fontsize=7.5, color="#0277bd", style="italic")

# 生物素（中心口袋）
ax_bio = fig.add_axes([0.290, 0.45, 0.085, 0.30])
draw_molecule(ax_bio, BIOTIN, atom_r=0.16)

# 3) 虚拟筛选流程（中）
ax_pipe = fig.add_axes([0.445, 0.10, 0.195, 0.80])
ax_pipe.set_xlim(0, 1); ax_pipe.set_ylim(0, 1); ax_pipe.axis("off")

steps = [
    ("3,311\napproved\ndrugs", "#ffe0b2", "#e65100"),
    ("virtual\nscreening\n(docking)", "#ffccbc", "#bf360c"),
    ("100 ns\nmolecular\ndynamics", "#ffab91", "#870000"),
]
yh, gap, top = 0.20, 0.14, 0.92
for i, (txt, fc, ec) in enumerate(steps):
    y = top - i * (yh + gap)
    ax_pipe.add_patch(FancyBboxPatch((0.1, y - yh), 0.8, yh,
                                     boxstyle="round,pad=0.02,rounding_size=0.03",
                                     fc=fc, ec=ec, lw=1.8))
    ax_pipe.text(0.5, y - yh / 2, txt, ha="center", va="center",
                 fontsize=8.5, fontweight="bold", color="#4e342e")
    if i < len(steps) - 1:
        ax_pipe.annotate("", xy=(0.5, y - yh - 0.012), xytext=(0.5, y - 0.018),
                         arrowprops=dict(arrowstyle="-|>", color="#616161", lw=2))

# 4) 右：ureido 对比（两个分子轴）+ 虚线
ax_bio2 = fig.add_axes([0.685, 0.34, 0.14, 0.42])
draw_molecule(ax_bio2, BIOTIN, atom_r=0.15)
ax_bio2.set_title("biotin\n(ureido ring)", fontsize=8, color="#2e7d32", pad=4)

ax_barb = fig.add_axes([0.855, 0.34, 0.14, 0.42])
draw_molecule(ax_barb, BARB, atom_r=0.15)
ax_barb.set_title("barbituric acid\n(malonylurea)", fontsize=8, color="#b71c1c", pad=4)

# 虚线等价线 + 标签（独立轴画在对比上方）
ax_comp = fig.add_axes([0.665, 0.28, 0.335, 0.06])
ax_comp.axis("off")
ax_comp.set_xlim(0, 1); ax_comp.set_ylim(0, 1)
ax_comp.annotate("", xy=(0.82, 0.4), xytext=(0.16, 0.4),
                 arrowprops=dict(arrowstyle="<->", linestyle="--",
                                 color="#616161", lw=1.8))
ax_comp.text(0.5, 0.75, "ureido-ring mimicry", ha="center", va="center",
             fontsize=9, color="#424242", style="italic", fontweight="bold")

# 机制轴（底部横向）
ax_mech = fig.add_axes([0.22, 0.015, 0.76, 0.05])
ax_mech.set_xlim(0, 1); ax_mech.set_ylim(0, 1); ax_mech.axis("off")
mech = ["SMVT ↑", "biotin uptake ↑", "ACC", "FASN",
        "de novo lipogenesis", "tumor growth ↑"]
n = len(mech)
for i, m in enumerate(mech):
    x0 = i / n + 0.005
    w = 1 / n - 0.02
    last = (i == n - 1)
    ax_mech.add_patch(FancyBboxPatch((x0, 0.15), w, 0.7,
                                     boxstyle="round,pad=0.01,rounding_size=0.02",
                                     fc="#f3e5f5" if last else "#eceff1",
                                     ec="#8e24aa" if last else "#607d8b", lw=1.2))
    ax_mech.text(x0 + w / 2, 0.5, m, ha="center", va="center", fontsize=7,
                 fontweight="bold" if last else "normal",
                 color="#4a148c" if last else "#37474f")
    if i < n - 1:
        ax_mech.annotate("", xy=(x0 + w + 0.018, 0.5), xytext=(x0 + w + 0.002, 0.5),
                         arrowprops=dict(arrowstyle="-|>", color="#90a4ae", lw=1.5))

# 保存
plt.savefig(OUT + ".png", dpi=100, facecolor="white")
plt.savefig(OUT + ".pdf", facecolor="white")
print(f"Saved: {OUT}.png / {OUT}.pdf")
print(f"PNG size: {Image.open(OUT + '.png').size}")
