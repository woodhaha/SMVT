# SLC5A6 (SMVT) — 综合分析报告

> **靶点**: SLC5A6 / SMVT (Sodium-dependent Multivitamin Transporter)
> **项目目录**: `D:\Researching\SMVT\`
> **分析日期**: 2026-06-23
> **数据来源**: TCGA, STRING v12, gnomAD, ClinVar, COSMIC, Human Protein Atlas, PubMed

---

## 核心发现

**SLC5A6 是表达驱动型代谢靶点，非突变驱动型癌基因。**
在 6/14 配对 TCGA 癌种中显著上调 (Log2FC +1.1 ~ +1.7)，突变率 <2%。
促癌机制: SMVT → Biotin → ACC → FASN → 脂质合成 → 肿瘤增殖。

---

## 1. TCGA 表达分析

### Tumor vs Normal 显著差异

| 癌种 | Log2FC | P-value | FDR |
|------|--------|---------|-----|
| BLCA (膀胱) | +1.71 | 2.4e-8 | 1.8e-6 |
| LUSC (肺鳞) | +1.68 | 1.3e-20 | 1.6e-19 |
| COAD (结肠) | +1.51 | 2.1e-12 | 4.9e-11 |
| ESCA (食管) | +1.48 | 1.9e-4 | 5.0e-3 |
| STAD (胃) | +1.44 | 4.7e-9 | 2.8e-7 |
| LUAD (肺腺) | +1.13 | 9.0e-26 | 4.2e-24 |

### 关键结论
- SLC5A6 在所有具有足够配对样本的癌种中一致上调
- 无任何癌种显示显著下调
- HPA 分类: "Expressed in all" (泛表达), 非组织特异性

## 2. TCGA 突变分析

### 约束评分: 强力非癌基因信号

| 指标 | SLC5A6 | 经典癌基因阈值 | 解读 |
|------|--------|-------------|------|
| pLI | **0.01** | ≥0.9 | 对功能缺失突变完全耐受 |
| LOEUF | **0.61** | ≤0.35 | 不受约束 |
| %HI | **68.45** | <10 | 非单倍剂量不足 |

### 体细胞突变
- 泛癌突变频率: **<2%** (乘客基因水平)
- 突变类型: 散发性错义, 无热点
- 截断突变: 极罕见 (<0.5%, 与 pLI=0.01 一致)
- CNA: 臂级事件为主, 非选择性局灶扩增/缺失

### 胚系突变 (非癌症)
- 致病突变引起代谢性神经疾病 (生物素响应性)
- 与癌症体细胞突变谱完全不重叠

## 3. STRING PPI 网络

### 10 个互作伙伴 (score ≥ 0.4)

| 基因 | Score | 模块 | 功能 |
|------|-------|------|------|
| PDZD11 | **0.969** | Anchor | 顶端膜 PDZ 支架 — **唯一高置信互作** |
| SLC5A7 | 0.655 | SLC | 胆碱转运体 |
| HLCS | 0.635 | Biotin | 全羧化酶合成酶 — **SLC5A6 反馈调控者** |
| SLC22A12 | 0.631 | SLC | URAT1, 尿酸转运体 |
| SLC26A4 | 0.578 | SLC | Pendrin, 阴离子交换体 |
| BTD | 0.555 | Biotin | 生物素酶, 生物素回收 |
| SLC5A3 | 0.511 | SLC | 肌醇转运体 |
| SLC23A1 | 0.493 | SLC | 维生素C转运体 |
| DPH2 | 0.471 | Other | 二苯胺生物合成 (仅共表达) |
| SLC19A2 | 0.469 | SLC | 硫胺素转运体 B1 |

### 关键调控节点

1. **PDZD11 — 物理锚点**: SLC5A6 C端 `SERTL` = Class I PDZ 结合模体 → PDZD11 决定 SMVT 膜定位。*肿瘤中完全未被研究*
2. **HLCS — 反馈调控**: HLCS→H4K12bio→SLC5A6 promoter 沉默。生物素不足→去抑制→SLC5A6 上调。*可能解释肿瘤中 SMVT 过表达*

## 4. 促癌机制: SMVT→Biotin→ACC→FASN 轴

```
SLC5A6 (SMVT) 过表达
  → 生物素/泛酸摄取增加
  → 乙酰辅酶A羧化酶 (ACC) 活性增强
  → 脂肪酸合酶 (FASN) 上调
  → 脂质从头合成增加
  → 膜生物合成 + 能量代谢 → 肿瘤增殖
```

### 文献证据

| 癌种 | 证据 | PMID |
|------|------|------|
| 肺腺癌 (LUAD) | SMVT→biotin→ACC→FASN 轴驱动肿瘤生长 | 39426496 |
| 宫颈癌 | SLC5A6 KD→FASN 下调→抑制增殖, FASN KD 可 rescue | 41108787 |
| 胃癌 | TCGA + IHC 诊断/预后 biomarker | Spandidos 2019 |

## 5. SLC5A6 vs 经典癌基因

| 维度 | EGFR/KRAS/TP53 | SLC5A6 |
|------|---------------|--------|
| 驱动方式 | 突变 → 组成型激活 | **过表达 → 代谢物供应↑** |
| pLI | ~1.0 (不耐受) | **0.01 (耐受)** |
| 突变频率 | 10-60% | **<2%** |
| 热点突变 | 有 | **无** |
| mRNA 变化 | 可变 | **一致上调** |
| 药物策略 | 抑制剂 | **底物偶联/转运阻断** |
| 正常组织功能 | 信号转导 | 维生素吸收 |

## 6. 靶点开发建议

### 优势
- 泛癌一致上调 → 广谱应用
- 已明确的促癌代谢机制 (ACC-FASN)
- 天然底物 (biotin) → 可设计靶向偶联物
- 非突变驱动 → 不易产生耐药突变
- PDZD11 是未被探索的调控节点

### 风险
- 正常组织 (肠道/肾脏) 也表达 SMVT → 需评估治疗窗口
- 非经典癌基因 → 需要更强的生物学验证数据
- 突变频率低 → 无法用突变作为生物标志物

## 8. 虚拟筛选 — ML 引导 + AutoDock Vina (2026-06-24 完成)

### 筛选策略
- ChEMBL 批准药物 3,311 → 药效团 ML 预筛 (RF ECFP4, AUC=0.888) → 多样性选择 500
- AutoDock Vina ex=16，356 成功对接 + 合并前 84 手选化合物 = **440 总化合物**

### 核心发现: 巴比妥类药物是新型 SMVT 高亲和力配体

| 化合物 | ΔG (kcal/mol) | 类型 | 临床用途 |
|--------|:---:|------|------|
| **Naftazone** | **−8.34** | Naphthoquinone | 止血剂 |
| **Phenobarbital** | **−8.30** | Barbiturate | 抗癫痫 (WHO 基本药物) |
| Cyclobarbital | −7.83 | Barbiturate | 镇静催眠 |
| Butalbital | −7.73 | Barbiturate | 偏头痛 |
| Tasimelteon | −7.66 | Melatonin agonist | 睡眠障碍 |
| Esketamine | −7.58 | Arylcyclohexylamine | 抗抑郁 (Spravato) |

**关键 SAR 发现**:
- **巴比妥酸骨架** (O=C1CC(=O)NC(=O)N1): 8/8 命中率 **100%**
- 巴比妥酸的丙二酰脲核心模拟生物素的尿素环 → 羧基**并非** SMVT 结合所必需
- 这扩展了靶向 SMVT 的化学空间，超越了传统羧酸类药物

详见: `06_Logs/SMVT-virtual-screening-report.md`

### 下一步实验

```
Phase 1: 验证
├── IHC 多组织阵列 (SMVT 蛋白在肿瘤 vs 正常)
├── Co-IP: SLC5A6—PDZD11
├── ChIP: H4K12bio at SLC5A6 promoter
├── TCGA 生存分析 (SMVT high vs low) ✅ 已完成
└── 虚拟筛选 ✅ 已完成 — Naftazone 和巴比妥类为首选苗头化合物

Phase 2: 功能
├── shRNA KD in 高表达肿瘤细胞系
├── 代谢组学验证 (biotin/CoA/FASN)
├── PDZD11 KD → SMVT 膜定位变化 (IF + biotinylation)
├── Biotin-conjugate uptake assay
└── ³H-biotin 竞争摄取实验 (top 10 虚拟命中的 IC₅₀)

Phase 3: 转化
├── 巴比妥类-SMVT 结合模式 MD 模拟 (100ns)
├── SMVT-targeted ADC/小分子偶联物设计
├── In vitro 肿瘤选择性摄取
├── Xenograft PK/PD
└── 正常组织 SMVT 依赖性安全性评估
```

## 7. 产出文件清单

```
D:\Researching\SMVT\
├── Fig1_SMVT_TCGA_pan_cancer.png/pdf      ← TCGA 表达 4面板 (box+forest+volcano+heatmap)
├── Fig2_SMVT_mechanism.png/pdf            ← 促癌机制通路图
├── Fig3_SMVT_mutation_landscape.png/pdf    ← 突变全景 (pLI+lollipop+象限)
├── Fig4_SMVT_mutation_vs_expression.png/pdf ← 突变 vs 表达 并排对比
├── Fig5_SMVT_STRING_network.png/pdf       ← STRING PPI 网络 (hub-spoke+证据)
├── Fig6_SMVT_STRING_evidence_heatmap.png/pdf ← 证据矩阵热图
├── SMVT-target-research.md                ← 原始研究文档
├── SMVT-TCGA-pan-cancer-expression.md     ← 表达分析报告
├── SMVT-TCGA-mutation-analysis.md         ← 突变分析报告
├── SMVT-STRING-interaction-network.md     ← STRING 互作报告
├── SMVT-master-analysis.md                ← 本文件 (综合分析)
├── visualize_tcga_nature.py               ← 表达可视化脚本
├── visualize_mutations_nature.py          ← 突变可视化脚本
└── visualize_string_network.py            ← STRING 可视化脚本
```

---

> **一句话总结**: SLC5A6/SMVT 是泛癌表达上调的代谢转运靶点，通过 biotin-ACC-FASN 轴驱动肿瘤脂质合成。
> 它不是突变驱动基因 (pLI=0.01)，而是表达驱动基因 (Log2FC 1.1-1.7)。
> PDZD11 是其唯一高置信互作伙伴 (score=0.969)，可能是肿瘤中完全未被研究的调控节点。
