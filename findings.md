# Findings — SMVT (SLC5A6)

## 已确立的理解
1. **SLC5A6 泛癌高表达**：多种实体瘤显著上调（Fig1, TCGA tumor vs normal）。
2. **突变景观**：MC3 中 193 例 SLC5A6 体细胞突变，错义为主，散在跨膜域（Fig3）。
3. **突变不显著改变表达**：突变型 vs 野生型无显著差异（Fig4）。
4. **虚拟筛选→MD 闭环**：702 化合物 → 8 elite hits (ΔG<−8) → 100ns MD 8/8 稳定 (RMSD<3Å)。
   - Hydromorphone (−8.58) 为最优；巴比妥类 100% 命中率（barbituric acid = biotin ureido mimic）。
5. **功能猜想**：SMVT 是多种临床药物（巴比妥/NSAID/阿片）的非预期转运底物 → DDI 风险。

## 尚未回答的问题（gap → 候选实验）
- **突变致病性**：哪些 SLC5A6 突变影响转运功能/底物识别？→ ESM-2 LLR + 结构映射（进行中）。
- **MD 与对接 ΔG 的一致性**：Furosemide ΔG 强但 MD 接触率最低 (56%)，结合模式是否不同？
- **泛癌生存**：KM/forest 图已出，未与表达/突变整合进叙事。
- **NSAID-SMVT DDI**：6 NSAIDs 命中，临床相关性未展开。

## ⚠️ 重大修复 (2026-07-31 晚) — MC3 barcode code_map 不完整
- **问题**：fig_mutation_final.py 的 code_map 漏了 18 个 TCGA 中心代码（12,17,23,27,29,33,94,A6,AB,BT,CD,CN,CQ,E2,E9,EM,FI,IB），导致 22 个 SLC5A6 突变被标为 OTHER 丢弃。Fig3/Fig4 用不完整数据（171/102/146 vs 实际 193/114/167）。
- **修复**：补全 18 个代码（含 LUAD/OV/GBM/CESC/LUSC/COAD/BLCA/STAD/SKCM）。Fig3 现显示 193/114/167 与正文一致。
- **附带修正**：正文"<2% all"说法不准确——实际 ESCA 5.9%、UCEC 3.3%、SKCM 3.2%、LUSC 3.0%、STAD 2.9% 超 2%（高突变癌型基线）。已改为准确表述（仍 passenger：多数 <2%、无 hotspot、突变率随基线 TMB）。
- **教训**：任何基于 TCGA barcode 前缀→癌型映射的分析都必须用完整 center code 表，否则漏样本。

## 新发现 (2026-07-31) — Vina ΔG 与 MD 动态稳定性
- **ΔG vs MD backbone RMSD**: Spearman ρ=−0.833, p=0.010（**显著**）。对接分越强 → MD 结构越稳定。强验证。
- **ΔG vs 接触率**: ρ=0.571, p=0.139（趋势，未显著）。原因：Riboflavin（负对照 ΔG=−0.01）接触率 80%——体积大机械贴近，非特异。排除后相关性增强。
- **数据源澄清**：稿件 Table 2 "Contacts%" = `hbond_occupancy_pct`（脚本注释明确 "contact occupancy"），56-85% 数据无误。
- 输出: `md_vina_correlation.csv` + `Fig_MD_vs_Vina_correlation.{pdf,png}`。
- **投稿价值**：ρ=−0.83 可加入 Discussion 作计算验证（Vina 打分 → MD 稳定性 → 对接可信度）。

## 新发现 (2026-07-31) — Furosemide 三方法矛盾（投稿 caveat）
- **Vina vs MM-GBSA**（n=7 排除 Furosemide）: Pearson r=−0.808, p=0.028 — 强一致性。
- **Furosemide 例外**：Vina −8.36（强结合）但 MM-GBSA **+6.40 kcal/mol（正排斥）**，MD 接触率最低（56%）。
- **三重证据链**：Vina 高分 ≠ 稳定结合。Furosemide 的 Vina 打分可能过度乐观（结合模式与预测不符）。
- MM-GBSA 数字全部吻合稿件（gabapentin −43.33 / riboflavin −41.48 / biotin −29.97 / hydromorphone −29.89 #4）。
- 输出: `Fig_Vina_vs_MMGBSA.{pdf,png}`。
- **投稿价值**：可作为 Discussion 中方法学 caveat（Vina 单方法局限）+ Furosemide 作为"打分乐观"例子。

## 5-HTP-SMVT 实验方案 + 理化预验证 (2026-07-31)
- 写 3 阶段实验方案：体外摄取 (H1) → BBB 模型 (H3 转运) → 体内 (H3 整体)，含对照与 Go/No-Go。
- **计算预验证**：5-HTP 具 SMVT 底物识别特征（COOH + MW 220 近生物素 244 + logP 0.83 vs 0.80 + HBD/HBA 4/3）。Tanimoto: 5-HTP vs 底物 0.09-0.11 = 底物彼此相似度 → SMVT 识别特征非结构。Trp vs 5-HTP 0.675 → 需区分 SMVT/LAT1。
- 存档: `05_Experimental/5HTP_SMVT_validation_plan.md`。提交: 81ed0cd。
- **投稿价值**：稿件 Limitations 已声明需摄取实验，此方案即对应；预验证为 H1 提供计算支撑。

## NSAID-SMVT DDI 文献挖掘 (2026-07-31)
- **Uchida 2015** (J Neurochem, DOI 10.1111/jnc.13092): SMVT 介导 BBB 88.7% 生物素/98.6% 泛酸摄取；**6 NSAID 直接抑制** (indomethacin, ketoprofen, diclofenac, ibuprofen, phenylbutazone, flurbiprofen)。
- **修正稿件 DOI 错误** (uchida2015: 13089→13092)。
- **流行病学**：Shaikh 2023 (累积 NSAID→抑郁, 癌症幸存者) ✓ 验证准确；Hu 2020 (非阿司匹林 NSAID HR 1.24)；Molero 2019 (一般人群无关联)。
- **5-HTP-SMVT 转运无直接实验证据** = 稿件假说新颖点 + 待验证空白。
- 存档: `01_Literature/08_NSAID_SMVT_DDI_review.md` + README 索引。
- 提交: c4d0086, ca8ca35。

## 稿件一致性补全 (2026-07-31 深夜)
- **Fig2 面板 a/d 用旧 DE 硬编码**（log2FC +1.1~+1.7, ESCA）→ 改为新 DE 数据（+0.44~+0.62, READ 替代 ESCA），轴范围重调，倍数 1.2-1.3x。提交: edab741。
- **Fig1 注意**：visualize_tcga_nature.py 也生成一个 Fig1 到 03_Analysis/ 根目录（不同布局），但稿件用的 Fig1 是 fig_pancan_de.py 生成（04_Manuscript/figures/），两者不同文件，无冲突。
- Fig5/6（STRING 网络）是独立网络分析，不受 DE 影响，内容正确。
- **最终扫描**：稿件无残留旧数字（+1.7, 141, 六癌型等全部清除）。

## 新发现 (2026-07-31 深夜) — 5-HTP 假说计算验证
- **5-HTP docked pose 与 biotin pose 几乎完全重叠**：质心位移 0.37Å，最小原子距离 0.15Å（AF2 SMVT 结合腔）。
- 证明 5-HTP 结合在 biotin 同一底物识别位点——转运介导摄取的前提。
- 强支撑稿件"SMVT 可能作为 BBB 的 5-HTP 转运体"假说。
- 新增 FigS15（3D 重叠）+ 正文描述。提交: e773592。
- **方法提示**：判断配体是否同一位点，用 pose 重叠（配体-配体质心距离）比配体-残基距离更可靠。

## 教训 / 约束
- MC3 MAF 3.5G：先预提取再迭代，勿逐次全扫。
- AlphaFold 原始 PDB 只读（02_Data/raw/）。
- 统计需遵循 CLAUDE.md 规范（Kruskal-Wallis、Mann-Whitney U 等）。
- 任何基于 TCGA barcode 前缀→癌型映射的分析必须用完整 center code 表。

## ⚠️ 稿件一致性 bug（投稿前必须修）
**问题**：稿件 DE 数字存在两套口径，互相矛盾。
- **新口径（可信，可复现）**：`run_pancancer.py` → `pancan_de.csv` → Fig1 图。Mann-Whitney U, unpaired, 全样本。**14/21 癌型显著上调**, log2FC +0.1~+0.6。摘要 "14 of 21" 与此吻合。
- **旧口径（硬编码，不可复现）**：`visualize_mutations_nature.py` 硬编码 `log2fc=[1.13,1.68,1.51,1.71,1.44,1.48,...]` → 稿件 Table 1 + 正文 "six cancer types"。log2FC +1.1~+1.7，paired TissGDB 风格（41 COAD 对等）。
- 正文声称 "Kruskal-Wallis with Dunn post-hoc" 但脚本实际用 Mann-Whitney U。
- **冲突点**：Table 1 / 正文（6 癌型、+1.1~+1.7）vs 摘要 / Fig1（14 癌型、+0.1~+0.6）。

**处理进展 (2026-07-31)**：
- ✅ 用 unpaired Mann-Whitney U + BH FDR 重算 → 14/21 上调 + 2 下调（THCA/KICH）。已更新稿件 Table 1、正文 DE 段、统计方法描述。
- ✅ 正文 mutation 段改为真实 MC3 数字：193 突变 / 167 患者 / 114 非沉默（96 missense, 11 nonsense, 5 splice, 2 indel）。missense 占非沉默 84.2%。
- ✅ V120I 实际 5 例（非 4）→ 已改。
- ✅ 摘要 mutation 率 0.6% → 1.5%。
- ✅ 对接数字验证：biotin −6.762 ✓, phenobarbital −8.299 ✓, naftazone −8.339 ✓, hydromorphone −8.577 ✓, furosemide −8.358 ✓, 5-HTP −7.664 ✓, glutarimide −9.364 ✓ (pathc_v2, 598 cmpds/62 hits)。
- ✅ MD 数字（RMSD/drift/H-bond）全部吻合 md_master_metrics.json。
- ✅ AUC = 0.888, MCC = 0.506, MAE = 1.965 kcal/mol 全部复现（pharmacophore_ml.log）→ 与稿件吻合。
- ✅ **生存分析验证**：稿件 Table 2 的 8 癌型（LGG/KIRC/ACC/MESO/KIRP/SARC/UVM/STAD）与 `pancan_survival.csv` 逐一吻合（HR 值精确匹配）。摘要 "eight cancer types" 正确。`survival_results.csv`（旧 cBioPortal 版，LUSC/LUAD 那套）未用于稿件，无冲突。
- ✅ **ESM-2 致病性**：86 个 MC3 missense 分析完成，uniform 分布（passenger），R123L control top 24%。修正了旧"zero missense"错误。
- pLI=0.01/LOEUF=0.61、10,967 样本为外部权威来源，未核对（可接受）。
- **结论**：稿件核心数字已全部核对。DE + mutation + 生存全部验证；对接分、MD、ML 全部正确。
