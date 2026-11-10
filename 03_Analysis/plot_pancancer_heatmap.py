#!/usr/bin/env python3
"""泛癌 SMVT (SLC5A6) 表达热图 — 图形摘要组件。

读 pancan_de.csv，取显著癌型（14 上调 + 2 下调），
画 Normal vs Tumor 双列热图（coolwarm 配色，红=高表达），
右侧标注 log2FC。输出 PNG(300dpi) + PDF(矢量)。
"""
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

CSV = r"D:\Researching\SMVT\03_Analysis\outputs\pancan_de.csv"
OUT = r"D:\Researching\SMVT\03_Analysis\outputs\pancancer_heatmap"

# ── 读取数据 ──
rows = []
with open(CSV, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if r["sig"].strip() == "*":
            rows.append({
                "cancer": r["cancer"],
                "n_mean": float(r["n_mean"]),
                "t_mean": float(r["t_mean"]),
                "log2FC": float(r["log2FC"]),
            })

# 按 log2FC 降序（上调在前，下调在后）
rows.sort(key=lambda r: -r["log2FC"])
labels = [r["cancer"] for r in rows]
n_mean = np.array([r["n_mean"] for r in rows])
t_mean = np.array([r["t_mean"] for r in rows])
fc = np.array([r["log2FC"] for r in rows])

M = np.column_stack([n_mean, t_mean])  # (N, 2): Normal | Tumor

# ── 出版级 rcParams ──
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.linewidth": 0.6,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "svg.fonttype": "none",
})

# ── 绘图 ──
fig, ax = plt.subplots(figsize=(3.2, 4.6))

# 颜色范围覆盖数据 + 少量 padding，中心取全局均值附近
vmin, vmax = 1.72, 2.72
vcenter = 2.15
norm = TwoSlopeNorm(vmin=vmin, vcenter=vcenter, vmax=vmax)

im = ax.imshow(M, cmap="coolwarm", norm=norm, aspect="auto")

ax.set_xticks([0, 1])
ax.set_xticklabels(["Normal", "Tumor"], fontsize=10, fontweight="bold")
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=8)

# 列分隔线（细）
ax.axvline(0.5, color="white", linewidth=1.2)

# 单元格内写表达值（log2 FPKM）
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center",
                fontsize=7, color="black",
                fontweight="bold" if j == 1 else "normal")

# 右侧标注 log2FC（红色上调 / 蓝色下调）
for i, v in enumerate(fc):
    c = "#c62828" if v > 0 else "#1565c0"
    ax.text(2.55, i, f"{v:+.2f}", ha="left", va="center",
            fontsize=7.5, color=c, fontweight="bold")

ax.set_xlim(-0.5, 3.1)

# 底部 log2FC 标题
ax.text(2.55, len(labels) - 0.15, "log2FC", ha="left", va="top",
        fontsize=7.5, color="gray", fontstyle="italic")

# 顶部标注：上调数量
n_up = int((fc > 0).sum())
ax.set_title(f"SMVT mRNA — {n_up}/{len(rows)} cancer types upregulated",
             fontsize=9.5, fontweight="bold", pad=8)

# colorbar（log2 表达）
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
cbar.set_label("log2 expression", fontsize=8)
cbar.ax.tick_params(labelsize=7)

# 去除多余脊线
for s in ("top", "right", "left", "bottom"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)

plt.tight_layout(pad=0.5)
plt.savefig(OUT + ".png", bbox_inches="tight", facecolor="white")
plt.savefig(OUT + ".pdf", bbox_inches="tight", facecolor="white")
print(f"Saved: {OUT}.png / {OUT}.pdf")
print(f"显著癌型 {len(rows)} 个：{n_up} 上调 + {len(rows)-n_up} 下调")
