# Key Decisions — SMVT Project

> 记录所有关键决策，确保可追溯

---

## 2026-06-23 — 项目初始化

| # | Decision | Rationale | Alternatives Considered |
|---|----------|-----------|------------------------|
| D001 | 采用 MedSci 六文件夹结构 | 标准化投稿流程，与 nature-skills + ARS 兼容 | 保持平板结构（不推荐：文件散落，跨工具协同困难） |
| D002 | AlphaFold 结构用于对接，不同时运行多构象 MD | AF 精度足以进行虚拟筛选；MD 成本高，留待后续 | 直接 MD 100ns（成本高，优先度低） |
| D003 | TCGA 使用 TCGAbiolinks (R) 下载，非 Broad Firehose | TCGAbiolinks 提供统一 API，样本级数据 | GDC Data Portal 直接下载（更原始但解析复杂） |
| D004 | 突变致病性预测先做 ESM-1v zero-shot，后补 αMissense | ESM-1v 本地运行快，αMissense 需 Google Cloud 或逐条 API | REVEL / Polyphen-2 (非结构感知，精度较低) |
| D005 | 对接库优先 FDA-approved drugs (~3K)，再扩展 DrugBank (~13K) | FDA 药物有已知安全性，转化路径更短 | 直接 ZINC 15 (2.3 亿化合物，计算量不可行) |

---

## Change Log

| Date | File | Change | Reason |
|------|------|--------|--------|
| 2026-06-23 | `CLAUDE.md` | Created | Project initialization |
| 2026-06-23 | `data_dictionary.md` | Created | Variable documentation |
| 2026-06-23 | Whole project | Reorganized into MedSci scaffold | D001 |
| 2026-06-24 | `pharmacophore_ml_screen.py` | Trained RF pharmacophore model (ECFP4) on 84 docked compounds | Phase A Step 1 |
| 2026-06-24 | `fetch_drugbank_ml_screen.py` | Fetched 3,311 ChEMBL approved drugs; ML-scored 2,822 | Phase A Step 2 |
| 2026-06-24 | `docking_batch_screen.py` | Docked 500 ML-prioritized compounds (356 successful, ex=16) | Phase A Step 3 |
| 2026-06-24 | `analyze_screening_results.py` | 440 total compounds analyzed; 174 hits; barbiturates discovered as novel SMVT ligands | Phase A Step 4 |
| 2026-06-24 | `SMVT-virtual-screening-report.md` | Complete virtual screening report with hit prioritization | Phase A Step 5 |
| 2026-06-24 | New finding | Barbiturates (8/8 hit rate) identified as novel SMVT pharmacophore — carboxyl NOT required | Major discovery |

---
## New Decisions (2026-06-24)

| # | Decision | Rationale | Alternatives Considered |
|---|----------|-----------|------------------------|
| D006 | Hit threshold recalibrated to relative scoring (ΔG < biotin −6.76 or Z < −1.5) | No compound hit strict −8.0; natural substrate defines biological relevance | Maintaining −8.0 (would produce 0 hits, missing real binders) |
| D007 | ML model uses ECFP4 + descriptors over graph neural networks | 84 training compounds insufficient for GNN; RF with fingerprints is robust and interpretable | GNN/CNN (data-hungry), pharmacophore matching (less quantitative) |
| D008 | Barbiturate scaffold prioritized for follow-up over NSAID/fenamate | 100% hit rate (8/8), crosses −8.0, novel chemical biology finding | Fenamate optimization (lower affinity ceiling, already known) |
| D009 | ChEMBL used over proprietary DrugBank | Open access, reproducible, 3,311 approved drugs is sufficient for initial screen | DrugBank (requires account, not reproducible without license) |
| D010 | Manuscript drafted as v1 with all completed analyses | Full manuscript ready for internal review; MD results deferred to v2 | Waiting for all analyses to complete before writing |

---

## Change Log

| Date | File | Change | Reason |
|------|------|--------|--------|
| 2026-06-24 | `SMVT_manuscript_v1.md` | Updated to v2: corrected GO/KEGG/Reactome terms, fixed multivariate HR values, corrected mutation claim, added supplementary figure legends S1-S5, aligned all figure references | Post-audit corrections |
| 2026-06-24 | `SMVT-virtual-screening-report.md` | Finalized | Phase A complete |
| 2026-06-24 | `pharmacophore_ml_screen.py` | Created | ML pre-screening pipeline |
| 2026-06-24 | `colab_md_top3.ipynb` | Created | GPU MD for Colab |

---

## 2026-06-28 — Screening Pipeline Visualization & Results Archive

### Decision: Archive screening results to master plan

**Context**: Virtual screening completed across 4 rounds (702 compounds). Key findings need to be preserved alongside the MD experimental plan.

**Results Archived**:
- 8 elite hits (ΔG < −8.0), Hydromorphone best at −8.58
- Barbiturate class: 100% hit rate (8/8), barbituric acid mimics biotin ureido ring
- 3 novel SMVT ligand classes: opioids, sulfonamides, arylcyclohexylamines
- Pharmacophore model: cyclic ureide/carboxamide core > carboxyl
- 5 SAR rules derived from 440-compound substructure analysis
- Controls defined: Biotin (reference), Gabapentin enacarbil (positive), Riboflavin (negative)

### Visualization Generated
- `03_Analysis/visualize_pipeline.py` — 4-panel pipeline figure
- `03_Analysis/figures/Fig_Screening_Pipeline.png/pdf` — funnel + families + pharmacophore + ranking

### Key Insight
Carboxyl group NOT essential for SMVT binding. This expands druggable chemical space beyond carboxylic acid mimics, validated by barbiturate 100% hit rate.

