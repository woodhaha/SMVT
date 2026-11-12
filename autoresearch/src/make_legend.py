"""Draw the color legend for the structure figures (matplotlib, run with base python)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIG = r"D:\Researching\SMVT\autoresearch\paper\figures"
RED = (0.714, 0.263, 0.259)
TEAL = (0.259, 0.580, 0.620)
GREEN = (0.180, 0.620, 0.267)

fig, ax = plt.subplots(figsize=(5.2, 1.5))
for i, (c, label) in enumerate([
    (RED, "COOH  —  SMVT anchor (free carboxyl)"),
    (TEAL, "CF3  —  affinity enhancer (x2)"),
    (GREEN, "7-position  —  conjugation handle"),
]):
    ax.add_patch(plt.Rectangle((0.02, 0.72 - i * 0.32), 0.14, 0.20, color=c))
    ax.text(0.20, 0.82 - i * 0.32, label, va="center", fontsize=11, family="sans-serif")
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
fig.savefig(FIG + "/structure_legend.png", dpi=200, bbox_inches="tight")
print("saved legend")
