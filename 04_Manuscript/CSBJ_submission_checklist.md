# CSBJ 投稿 Checklist（Editorial Manager）

> 状态：✅ 已完成 · ⚠️ 待处理 · ❌ 缺失

## A. 核心投稿文件

- [x] **稿件 PDF** — `SMVT_manuscript.pdf`（42 页，编译 0 错误）
- [x] **Cover letter** — `cover_letter_CSBJ.md`（计算化学 + 药物发现 framing）
- [ ] **Title page（单独文件）** — 作者名 / 单位 / ORCID / 通讯邮箱。稿件里作者信息在 `\author`/`\affil`，但 Elsevier 要单独 title page。
- [x] **Abstract + Keywords** — 已在稿件内（非结构化摘要 + 10 关键词）
- [ ] **Highlights（3–5 条要点）** — Elsevier 投稿系统强制要求，每条约 85 字符。
- [ ] **Graphical abstract** — 可选但强烈建议（提高 editor 初筛通过率）。
- [x] **主图 + 补充图** — 22/24 已统一（见下）。
- [x] **表格** — Table 1–4 在稿件内。
- [x] **补充材料** — FigS1–S16 + 图注。

## B. 声明项（稿件 Declarations 段已含）

- [x] Ethics approval（公开去标识数据，无需伦理）
- [x] Funding（无资助）
- [x] Conflict of interest（无）
- [x] Author contributions（CRediT 式：Z.H. 构思/分析/撰写）
- [x] Data availability（TCGA/GDC、GEO、gnomAD、ChEMBL、AlphaFold、RCSB）
- [x] Code availability（Zenodo/Figshare upon acceptance + 审稿期可提供）

## C. Editorial Manager 系统项

- [ ] **通讯作者 ORCID** — 确认用 `0009-0005-1818-4789`（活跃，已绑 Springer；Elsevier 需另授权）
- [ ] **Suggested reviewers（3–5 名）** — 需提供姓名/机构/邮箱，建议选自 SLC 转运体、biotin 代谢、计算药物发现领域。
- [x] **Manuscript type** — Research Article
- [x] **Keywords（系统用）** — 与稿件一致
- [ ] **Funding source（系统勾选）** — 填 "None"

## D. 格式 / 一致性

- [ ] **LaTeX 模板** — ⚠️ 稿件用 Springer `sn-jnl` class，CSBJ 是 Elsevier。初投可接受 PDF，但**接收后要重排成 Elsevier 格式**（或改投前用 Elsevier LaTeX/Word 模板重排更稳妥）。
- [x] 图分辨率 ≥300 dpi
- [x] 参考文献 — 24 条，全部 CrossRef/PMID 验证过（含 tripathi2025）
- [ ] **图正文引用核对** — 确认 FigS1–S16 每个都在正文被引用、无孤儿。

## E. 待补的 2 张图（投稿前必须）

- [ ] **FigS9（结合模式）** — 需 PyMOL 渲染（`MM-GBSA Analysis/_viz_binding_mode.py`）
- [ ] **FigS15（5-HTP pose 重叠 + RMSD）** — 需 docking pose 数据（`05_PDB/validation/pose_comparison.py`）

## 建议 Highlights 草稿（供直接套用）

- SLC5A6 (SMVT) is upregulated in 14/21 cancer types yet carries only passenger mutations
- Survival associations (HR 1.8–4.0) survive multivariable adjustment for age, sex, and stage
- Virtual screening predicts barbiturates and glutarimides as SMVT ligands (100% barbiturate hit rate)
- 5-HTP and NSAIDs dock in the biotin site, motivating a testable SMVT–serotonin axis
- Cross-structure re-docking against cryo-EM benchmarks the AF2-based screening pipeline
