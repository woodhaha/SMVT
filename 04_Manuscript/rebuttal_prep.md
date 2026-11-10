# Rebuttal 预演 — SMVT/CBM（返修备用）

> 用途：投稿后若进 major revision，审稿人最可能问的三问 + 反驳逻辑。
> 原则：每条反驳均用稿件现有文字支撑，**不引入新数据/新实验**（无湿实验条件）。
> 定位：纯计算预测 + 可检验假设（CBM 主场）。

---

## Q1 — 巴比妥/戊二酰亚胺缺羧基，不能被转运，top hits 无意义

**反驳要点**
- 我们从没声称巴比妥/戊二酰亚胺**可转运**，全文明确标为 `pharmacophore scaffold`（predicted）。
- 引 tripathi2025 + zhen2026 说明羧基是转运必需（稿件第 307/391 行已写）。
- 真正的"预测底物"是**含羧基**的 5-HTP 与 NSAID，且复现了实验报道的 SMVT 抑制（uchida2015：六个 NSAID 抑制 biotin/泛酸摄取）。
- 巴比妥的价值是**药化起点**：malonylurea 模拟 biotin 的 ureido 环（已验证的氢键供体/受体模体），非转运声明。

**稿件对应位置**：摘要、Results 第 307 行、Discussion 第 389/391 行、结论。

---

## Q2 — MD 没 Na⁺ 没膜，SMVT 是 Na⁺ 偶联转运体，稳定性结果不可靠

**反驳要点**
- 已如实披露（第 397 行，并已补 zhen2026 对照）。
- 关键事实：**cryo-EM 结构本身也没解析出 Na⁺ 密度**——zhen2026 是同源放置 2 个 Na⁺（Na2/Na3，靠 SLC5 家族保守位点）。
- 我们 re-docking 用的是 cryo-EM **occluded 态 = Na⁺ 稳定构象**，羧基配位几何已部分涵盖。
- MD 用途是**固定口袋内的稳定性排序**，不算转运能量学；Na⁺ 主要影响羧基配位几何，由 cryo-EM 口袋兜底。
- Na⁺/膜 MD 已承诺为后续工作（limitation 内）。

**稿件对应位置**：Discussion 第 397 行（limitation 第四点）。

---

## Q3 — 标准套路（TCGA + Vina + OpenMM），计算新颖性在哪

**反驳要点**
- 贡献不是单个方法，是**整合 + 基准**：
  - (a) SLC5A6 首个泛癌**表达驱动型**脆弱点刻画（非单纯筛选练习）；
  - (b) **cross-structure re-docking 基准**——F1=91% 口袋重叠、hydromorphone 构象依赖性，是对 AF2 筛选的**可迁移方法学警示**；
  - (c) ML 药效团筛选（AUC 0.888）与 docking 整合。
- 可迁移结论：AF2 口袋对特定骨架（吗啡喃）会过拟合，跨 cryo-EM 构象 re-docking 应成为 AF2 筛选的标准步骤。

**稿件对应位置**：摘要、Results 第 343 行（re-docking）、第 281 行（ML 筛选）。

---

## Q4 — SMVT 作为药物靶点已有大量研究（biotin 偶联物递送），这篇有什么新？

**反驳要点**
- 前人工作是**递送**：用 SMVT 把 biotin-drug 偶联物带进细胞（2015 综述 PMC4406285；gabapentin enacarbil 递送前药）。
- 我们做的是**抑制/重定位**：筛选已批准药结合 SMVT 维生素位点，方向相反。
- 且用刚解出的 cryo-EM（zhen2026，2026-06）做跨构象 re-docking 基准——递送文献用同源模型做不到。
- 递送 vs 抑制已在 discussion 显式区分（新增句）。

**稿件对应位置**：Discussion "First" 段（递送 vs 抑制区分句）、摘要结论。

---

## 通用退让话术（审稿人要求湿实验时）

> "We fully agree that experimental validation is the definitive test of these predictions. In the absence of the facility for radiolabeled uptake or biophysical binding assays, we have (i) constrained all ligand claims as *predicted*, (ii) anchored them to the experimentally validated carboxyl-recognition chemistry (tripathi2025, uchida2015, zhen2026), and (iii) provided a concrete, ranked list of competition-uptake and SPR/ITC experiments (Supplementary/Future directions) for direct follow-up."

---

*注：本文件是预演草稿，返修时按 EM/期刊格式改写，并将稿件行号替换为终稿行号。*
