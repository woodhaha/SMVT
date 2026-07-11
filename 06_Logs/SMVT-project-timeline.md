# SMVT 项目时间线 — 最早步骤

> 项目代号: SMVT · 基因: SLC5A6 · 蛋白: Na⁺-依赖多维生素转运体
> 记录从零启动到虚拟筛选完成的完整研究步骤

---

## 2026-06-23 — 项目启动日

### Step 1: 靶点选择与文献调研
- 文件: `06_Logs/SMVT-target-research.md`
- 选择 SLC5A6/SMVT 作为靶点
- 唯一上市 SMVT 靶向药物: Gabapentin enacarbil (Horizant®)
- 核心发现: SMVT 在多种肿瘤中过表达（乳腺癌、前列腺癌、胃癌、胶质瘤）
- NSAIDs 被识别为 SMVT 抑制剂
- 2024-2025 活跃管线: Biotin 偶联策略（小分子 + 纳米粒 + 放射治疗）
- 工具: DrugBank-Database Skill + WebSearch

### Step 2: 结构准备
- 文件: `06_Logs/SMVT-structure-preparation-report.md` + `minimization.log`
- PDB 实验结构不存在（膜蛋白，结晶困难）→ AlphaFold v6 (AF-Q9Y289-F1)
- pLDDT=79.4, 635 残基, 4,824 原子
- 流水线: 下载 → Biopython 去杂 → OpenMM 加氢(pH 7.0) → Amber14SB+GBSA 能量最小化
- 最终能量 −55,275 kJ/mol (降幅 15,405 kJ/mol)
- 输出: `02_Data/cleaned/AF-Q9Y289-F1_prepared.pdb`
- 工具: AlphaFold DB, Biopython, OpenBabel, OpenMM

### Step 3: TCGA 泛癌表达分析 (Fig1)
- 文件: `06_Logs/SMVT-TCGA-pan-cancer-expression.md`
- SLC5A6 在 6/14 配对癌种中一致上调（Log2FC +1.1 ~ +1.7）
- 显著癌种: BLCA(+1.71), LUSC(+1.68), COAD(+1.51), ESCA(+1.48), STAD(+1.44), LUAD(+1.13)
- 无任何癌种下调 → **表达驱动型**靶点
- 工具: TCGAbiolinks (R), visualize_tcga_nature.py

### Step 4: TCGA 突变分析 (Fig3-4)
- 文件: `06_Logs/SMVT-TCGA-mutation-analysis.md`
- 突变率 <2% (乘客基因水平), pLI=0.01, LOEUF=0.61
- 无热点突变，截断突变极罕见
- 胚系致病突变 → 代谢性神经疾病（生物素响应性），与癌症体细胞突变谱不重叠
- 工具: gnomAD, ClinVar, COSMIC, visualize_mutations_nature.py

### Step 5: STRING PPI 网络 (Fig5-6)
- 文件: `06_Logs/SMVT-STRING-interaction-network.md`
- PDZD11 唯一高置信互作 (score=0.969) — 顶端膜 PDZ 支架，**肿瘤中完全未被研究**
- HLCS 反馈调控: H4K12bio → SLC5A6 启动子沉默 → 可能解释肿瘤中过表达
- SLC 家族互作: SLC5A7(胆碱), SLC22A12(尿酸), SLC26A4(阴离子), SLC5A3(肌醇), SLC23A1(维C), SLC19A2(硫胺素)
- 工具: STRING v12, visualize_string_network.py

### Step 6: 综合分析 — 核心假说
- 文件: `06_Logs/SMVT-master-analysis.md`
- **核心假说**: SLC5A6 过表达 → Biotin 摄取↑ → ACC → FASN → 脂质合成 → 肿瘤增殖
- 文献证据: PMID 39426496 (LUAD), PMID 41108787 (宫颈癌)
- **关键诊断**: SMVT 是**表达驱动型**代谢靶点，非突变驱动型癌基因

### Step 7: 关键决策
- 文件: `06_Logs/decisions.md`
- D001: 采用 MedSci 六文件夹标准化结构
- D002: AlphaFold 结构做对接，不先跑 MD（成本高，优先级低）
- D003: TCGA 用 TCGAbiolinks 下载
- D004: 突变致病性预测先 ESM-1v 后 αMissense
- D005: 对接库优先 FDA 获批药物 (~3K)，转译路径短

---

## 2026-06-24 — 虚拟筛选阶段

### Step 8: ML 药效团预筛
- 文件: `03_Analysis/pharmacophore_ml_screen.py`
- 84 个手选化合物 → RF (ECFP4 + 分子描述符) → AUC=0.888
- 预筛 ChEMBL 3,311 获批药物 → 2,822 个可评分化合物
- 决策 D007: RF+ECFP4 而非 GNN（84 样本不够深度学习）

### Step 9: 虚拟筛选执行
- 文件: `03_Analysis/docking_batch_screen.py`, `dock_parallel.py`
- 500 ML 优先化合物 → AutoDock Vina (ex=16)
- 356 成功对接 + 84 手选 = **440 总化合物**
- 多轮补充: dock_leftover_fda.py, dock_leftover_v2.py, docking_round3.py

### Step 10: 命中分析与 SAR
- 文件: `03_Analysis/analyze_screening_results.py`, `03_Analysis/sar_deep_analysis.py`
- **8 个 Elite Hits** (ΔG < −8.0): Hydromorphone(−8.58), Furosemide(−8.36), Naftazone(−8.34), Phenobarbital(−8.30)...
- **巴比妥类 100% 命中率** (8/8): 巴比妥酸骨架模拟生物素尿素环
- **关键发现**: 羧基并非 SMVT 结合所必需 → 扩展了靶向化学空间
- 其他命中类别: NSAID(6/12=50%), Opioid(2/5=40%), Sulfonamide(1/6=17%)

### Step 11: 虚拟筛选报告
- 文件: `06_Logs/SMVT-virtual-screening-report.md`
- 完整 SAR 规则 (5 条)
- 对照验证: Biotin(−6.76, 底物), Gabapentin enacarbil(−6.63, 阳性), Riboflavin(−0.01, 阴性)
- 决策 D008: 巴比妥类优先跟进
- 决策 D010: 稿件 v1 完成，MD 结果留到 v2

---

## 2026-06-28 — MD 实验准备 + TCGA 数据复验

### Step 12: TCGA 表达数据交叉验证
- 文件: `06_Logs/SMVT-TCGA-expression-validation.md`
- TissGDB 二次查询 (6/23 + 6/28) — 6 种癌种数据完全一致
- 本地 CRCproject TCGA 数据 (520 样本) — SLC5A6 检出
- HPA + 3 篇 PubMed 文献三线交叉验证
- ESCA 信号弱 (FDR=0.005, 仅 10 对配对样本)
- 8 种额外癌种无显著差异 (KIRC/KIRP/LIHC/BRCA/THCA/PRAD/HNSC/KICH)
- 可视化: `Fig_SMVT_TCGA_Validated.{png,pdf}` (4 面板: box+forest+volcano+validation)

### Step 13: MD 协议开发
- Esketamine 试点调试
- **关键发现**: 约束必须先加再最小化（v3/v5），不能先最小化再加约束（v7/v8 bug）
- Vina 对接姿势有 3 处严重空间冲突 (2.08Å=61% vdW)
- 4 阶段分级最小化验证: PE 从 7.6T 降至 −885K
- 执行策略: Colab Free T4, 1 compound/day

---

## 项目文件索引

| 阶段 | 关键文件 |
|------|---------|
| 靶点研究 | `06_Logs/SMVT-target-research.md` |
| 结构准备 | `06_Logs/SMVT-structure-preparation-report.md` |
| 表达分析 | `06_Logs/SMVT-TCGA-pan-cancer-expression.md` |
| 突变分析 | `06_Logs/SMVT-TCGA-mutation-analysis.md` |
| PPI 网络 | `06_Logs/SMVT-STRING-interaction-network.md` |
| 综合分析 | `06_Logs/SMVT-master-analysis.md` |
| 虚拟筛选 | `06_Logs/SMVT-virtual-screening-report.md` |
| 决策记录 | `06_Logs/decisions.md` |
| 主计划 | `Plans.md` |
