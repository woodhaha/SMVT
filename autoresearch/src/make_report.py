"""Generate progress report (bar chart + HTML) for SMVT drug-design autoresearch."""
import os, csv, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "to_human")

# collect all results
rows = []
for csvf in glob.glob(os.path.join(BASE, "experiments", "*", "*_results.csv")):
    with open(csvf) as f:
        for r in csv.DictReader(f):
            try:
                aff = float(r["best_affinity"])
            except ValueError:
                continue
            if aff < 0 and r.get("ok") == "True":
                rows.append((r["name"], aff))

# dedupe by name (keep most negative)
best = {}
for name, aff in rows:
    if name not in best or aff < best[name]:
        best[name] = aff

top = sorted(best.items(), key=lambda x: x[1])[:20]

# bar chart
plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "DejaVu Sans"],
                     "font.size": 9})
names = [n for n, _ in top][::-1]
vals = [v for _, v in top][::-1]
colors = ["#B64342" if v < -9.0 else "#0F4D92" if v < -8.5 else "#767676" for v in vals]
fig, ax = plt.subplots(figsize=(8, 6))
ax.barh(names, vals, color=colors)
ax.axvline(-8.84, color="#42949E", ls="--", lw=1.2, label="prior best transportable (−8.84)")
ax.axvline(-9.36, color="#9A4D8E", ls="--", lw=1.2, label="glutarimide (−9.36, carboxyl-free)")
ax.set_xlabel("Vina score (kcal/mol)")
ax.set_title("SMVT carboxyl-bearing ligands — top 20 by docking score")
ax.legend(frameon=False, fontsize=8)
fig.tight_layout()
png = os.path.join(OUT, "top20_scores.png")
fig.savefig(png, dpi=200)
plt.close(fig)

# HTML report
headline = top[0]
html = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>SMVT Drug Design — Autoresearch Report</title>
<style>body{{font-family:Arial,sans-serif;max-width:900px;margin:2rem auto;padding:0 1rem;line-height:1.5}}
h1{{color:#0F4D92}} h2{{color:#333;border-bottom:1px solid #ddd;padding-bottom:.3rem}}
table{{border-collapse:collapse;width:100%;font-size:14px}}
th,td{{border:1px solid #ddd;padding:6px 10px;text-align:left}}
th{{background:#f2f2f2}} .hit{{background:#fff3f3;font-weight:bold}}
.badge{{display:inline-block;background:#B64342;color:#fff;padding:2px 8px;border-radius:4px;font-size:13px}}
</style></head><body>
<h1>SMVT 药物设计 — 自主研究进度报告</h1>
<p><b>研究问题</b>：能否找到 <b>含自由羧基</b>（可转运）的 SMVT 配体，Vina 分数 < −8.84（超越当前最佳）？</p>
<p><b>答案</b>：<span class="badge">是</span> — <b>{headline[0]}</b> = <b>{headline[1]:.2f} kcal/mol</b>，3 次独立 Vina 种子复现（−9.68~−9.70），打破 −8.84 天花板，也超过无羧基的戊二酰亚胺 −9.36。</p>
<img src="top20_scores.png" style="width:100%;max-width:800px">
<h2>完整 SAR 轨迹（5 实验 ~75 化合物）</h2>
<table><tr><th>实验</th><th>化学型</th><th>最佳</th><th>结论</th></tr>
<tr><td>E1</td><td>芳香氨基酸</td><td>1-萘丙氨酸 −7.74</td><td>次级线索，封顶 ~−7.7</td></tr>
<tr><td>E2</td><td>新笼</td><td>金刚烷二羧酸 −8.48</td><td>替代笼全劣</td></tr>
<tr><td>E3</td><td>羧基+小基团</td><td>COOH-CF3 −9.02</td><td>首破 −8.84</td></tr>
<tr><td>E4</td><td>CF3 变体</td><td>bisCF3 −9.69</td><td>突破</td></tr>
<tr><td>E5</td><td>确认+推进</td><td>bisCF3 −9.68（3 seeds）</td><td>确认最优</td></tr>
</table>
<h2>核心发现</h2>
<p><b>三氟甲基（CF3）是金刚烷笼上的强亲和增强剂</b>——同时是吸电子（强化羧基锚定）+ 疏水（填口袋）+ 小（无位阻）。CF3 数量最优：1 个(−9.04) &lt; <b>2 个(−9.68)</b> &gt; 3 个(−8.03)。</p>
<h2>下一步（验证，非更多 docking）</h2>
<ol><li>cryo-EM 26va 交叉对接验证 COOH-bisCF3</li><li>SMVT–bisCF3 的 MM-PBSA / 膜 MD</li><li>湿实验：biotin 竞争摄取实验</li></ol>
</body></html>"""
html_path = os.path.join(OUT, "report.html")
with open(html_path, "w", encoding="utf-8") as f:
    f.write(html)
print("saved:", png)
print("saved:", html_path)
print(f"\ntop 5:")
for n, v in top[:5]:
    print(f"  {n:30s} {v:8.2f}")
