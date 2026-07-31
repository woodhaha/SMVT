# Findings — SMVT (SLC5A6)

## 已确立的理解
1. **SLC5A6 泛癌高表达**：多种实体瘤显著上调（Fig1, TCGA tumor vs normal）。
2. **突变景观**：MC3 中 193 例 SLC5A6 体细胞突变，错义为主，散在跨膜域（Fig3）。
3. **突变不显著改变表达**：突变型 vs 野生型无显著差异（Fig4）。
4. **虚拟筛选→MD 闭环**：702 化合物 → 8 elite hits (ΔG<−8) → 100ns MD 8/8 稳定 (RMSD<3Å)。
   - Hydromorphone (−8.58) 为最优；巴比妥类 100% 命中率（barbituric acid = biotin ureido mimic）。
5. **功能猜想**：SMVT 是多种临床药物（巴比妥/NSAID/阿片）的非预期转运底物 → DDI 风险。

## 尚未回答的问题（gap → 候选实验）
- **突变致病性**：哪些 SLC5A6 突变影响转运功能/底物识别？→ αMissense/ESM-1v + 结构映射。
- **MD 与对接 ΔG 的一致性**：Furosemide ΔG 强但 MD 接触率最低 (56%)，结合模式是否不同？
- **泛癌生存**：KM/forest 图已出，未与表达/突变整合进叙事。
- **NSAID-SMVT DDI**：6 NSAIDs 命中，临床相关性未展开。

## 教训 / 约束
- MC3 MAF 3.5G：先预提取再迭代，勿逐次全扫。
- AlphaFold 原始 PDB 只读（02_Data/raw/）。
- 统计需遵循 CLAUDE.md 规范（Kruskal-Wallis、Mann-Whitney U 等）。
