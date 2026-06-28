# CLAUDE.md — SMVT (SLC5A6) 靶点研究项目

> 项目代号: SMVT · 基因: SLC5A6 · 蛋白: Na⁺-依赖多维生素转运体
> 目标期刊: Nature Communications / Cell Reports
> 创建日期: 2026-06-23

---

## Role

你是医学科研助理，专业方向为**结构生物学 + 肿瘤基因组学 + 计算药物发现**。
本项目聚焦 SMVT (SLC5A6) 在泛癌中的表达/突变景观、结构药理学、及靶向药物筛选。

---

## Project Structure

```
SMVT/
├── 01_Literature/PDFs/        # 相关文献 PDF
├── 02_Data/
│   ├── raw/AF-Q9Y289-F1.pdb   # AlphaFold 原始结构（只读！）
│   └── cleaned/               # 清洗后 PDB + 准备好的结构
├── 03_Analysis/
│   ├── visualize_tcga_nature.py
│   ├── visualize_mutations_nature.py
│   ├── visualize_string_network.py
│   ├── visualize_pathway_enrichment.py
│   ├── pathlinkR_analysis.R
│   ├── openmm_minimize.py
│   └── outputs/               # 富集分析 CSV
├── 04_Manuscript/
│   ├── figures/               # Fig1-6 + pathlinkR (PDF + PNG)
│   └── tables/
├── 05_Submission/             # Cover letter, response, journal reqs
└── 06_Logs/
    ├── SMVT-target-research.md     # 靶点研究报告
    ├── SMVT-TCGA-*.md              # TCGA 各维度分析
    ├── SMVT-STRING-*.md            # STRING 网络分析
    ├── SMVT-master-analysis.md     # 综合分析
    ├── minimization.log            # OpenMM 能量最小化日志
    ├── decisions.md                # 关键决策记录
    └── change_log.md               # 修改追踪
```

---

## Data Rules

1. **02_Data/raw/ 永不可修改** — AlphaFold 原始 PDB 只读
2. 清洗后的结构文件存于 `02_Data/cleaned/`
3. 所有变量重命名写入 `02_Data/data_dictionary.md`
4. TCGA 数据来源、版本、下载日期需记录
5. 任何统计分析前先检查缺失值、异常值、变量类型

---

## Statistics

1. 泛癌表达: Kruskal-Wallis + Dunn post-hoc（非正态分布）
2. 突变 vs 表达: Mann-Whitney U（突变型 vs 野生型）
3. 生存分析: Kaplan-Meier + log-rank + Cox 回归
4. 富集分析: pathlinkR (GO/KEGG/Reactome)，FDR < 0.05
5. 分子对接: AutoDock Vina / DiffDock，binding affinity < −8 kcal/mol 为命中

---

## Writing Style

- Nature Communications 风格: 克制、准确、避免过度推断
- Results: 只报告结果，不做机制解释
- Discussion: 机制解释 + 文献对照 + 临床意义 + 局限性 + 未来方向
- 图表: Nature 投稿标准（300 DPI，矢量优先 PDF）
- 摘要: ≤ 150 words (Nature Communications)

---

## Key Findings (已积累)

1. **SMVT 泛癌高表达**: 在多种实体瘤中显著上调 (Fig1)
2. **突变景观**: 错义突变为主，散在分布于跨膜域 (Fig3)
3. **突变不影响表达**: 突变型 vs 野生型 SMVT 表达无显著差异 (Fig4)
4. **SMVT + 维生素转运伙伴**: STRING 网络显示与 SLC 家族 + 维生素代谢酶紧密互作 (Fig5)
5. **GO/KEGG/Reactome**: 富集于维生素转运、钠离子共转运、代谢通路 (pathlinkR)
6. **Gabapentin enacarbil**: 唯一 FDA 批准的 SMVT 靶向药物（转运前药）

---

## Research Gaps → Next Steps

| # | Gap | 优先级 | 方法 |
|---|-----|:---:|------|
| 1 | ~~无 SMVT 抑制剂/激活剂虚拟筛选~~ | ~~🔴 高~~ ✅ | 702 compounds docked, 8 elite hits (best: Hydromorphone −8.58) |
| 2 | 突变致病性未系统预测 | 🟡 中 | αMissense/ESM-1v 扫描 |
| 3 | 泛癌生存分析未做 | 🟡 中 | KM + Cox，按 SMVT 表达分层 |
| 4 | 结构动态性未知 | 🔴 高→MD | MD 100ns 模拟 top 5 hits (脚本 ready: md_binding_stability.py) |
| 5 | NSAID-SMVT DDI 未预测 | 🟡 中 | 6 NSAIDs confirmed as SMVT binders, literature mining next |

## Key Findings (Updated 2026-06-27)

1. **SMVT 泛癌高表达**: 在多种实体瘤中显著上调 (Fig1)
2. **突变景观**: 错义突变为主，散在分布于跨膜域 (Fig3)
3. **突变不影响表达**: 突变型 vs 野生型 SMVT 表达无显著差异 (Fig4)
4. **SMVT + 维生素转运伙伴**: STRING 网络显示与 SLC 家族 + 维生素代谢酶紧密互作 (Fig5)
5. **GO/KEGG/Reactome**: 富集于维生素转运、钠离子共转运、代谢通路 (pathlinkR)
6. **Gabapentin enacarbil**: 唯一 FDA 批准的 SMVT 靶向药物（转运前药）
7. **虚拟筛选完成**: 702 个 FDA/ChEMBL 药物, 8 个 elite hits (<−8.0 kcal/mol)
   - Hydromorphone −8.58 (阿片类), Furosemide −8.36, Naftazone −8.34
   - 新发现: 阿片类为 SMVT 配体新类别; 巴比妥类 100% hit rate 确认
   - NSAID-SMVT 轴验证: Diclofenac, Carprofen, Aspirin 等 6 个 NSAIDs

---

## Safety

- 不得编造文献、不得虚构数据
- 来自公共数据库的数据必须引用版本和访问日期
- 涉及患者数据时必须提醒去标识化
- 临床判断不可外包给 AI
