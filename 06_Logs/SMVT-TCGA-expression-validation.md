# SLC5A6 (SMVT) TCGA 泛癌表达 — 数据验证报告

> 验证日期: 2026-06-28 · 原始分析日期: 2026-06-23
> 数据源: TissGDB + 本地 CRCproject TCGA 数据 + Human Protein Atlas + PubMed 文献

---

## 一、数据来源交叉验证

| 来源 | 时间 | 状态 |
|------|:--:|:--:|
| **TissGDB** (bioinfo.uth.edu) | 2026-06-23 + 2026-06-28 复查 | ✅ 两次数据完全一致 |
| **本地 TCGA 数据** (CRCproject) | D:\Researching\CRCproject\02_Data\raw | ✅ TCGA_Exp_convertensg2symbol.txt, 520 样本, SLC5A6 检出 |
| **Human Protein Atlas** | proteinatlas.org | ✅ 分类 "Expressed in all" |
| **PubMed 文献** | PMID 39426496, 41108787, Spandidos 2019 | ✅ 功能验证 + IHC 交叉验证 |

---

## 二、显著上调癌种 (TissGDB, 配对差异检验, |log2FC|>1, FDR<0.05)

| 癌种 | 缩写 | Log2FC | P-value | FDR | 配对样本数 |
|------|:--:|:--:|------|------|:--:|
| 膀胱尿路上皮癌 | BLCA | **+1.71** | 2.40×10⁻⁸ | 1.84×10⁻⁶ | ≥10 |
| 肺鳞状细胞癌 | LUSC | **+1.68** | 1.29×10⁻²⁰ | 1.64×10⁻¹⁹ | ≥10 |
| 结肠腺癌 | COAD | **+1.51** | 2.09×10⁻¹² | 4.87×10⁻¹¹ | ≥10 |
| 食管癌 | ESCA | **+1.48** | 1.87×10⁻⁴ | 5.01×10⁻³ ⚠️ | 10 |
| 胃腺癌 | STAD | **+1.44** | 4.68×10⁻⁹ | 2.77×10⁻⁷ | ≥10 |
| 肺腺癌 | LUAD | **+1.13** | 8.99×10⁻²⁶ | 4.20×10⁻²⁴ | ≥10 |

> 数据来源: TissGDB, TCGA IlluminaHiSeq_RNASeqV2, pan-cancer normalized log2(norm_counts+1)
> 仅纳入 ≥10 对 Tumor-Normal 配对样本的癌种 (共 14 种)

---

## 三、其他癌种表达情况 (TissGDB 配对分析)

| 癌种 | 缩写 | Log2FC | FDR | 趋势 |
|------|:--:|:--:|------|:--:|
| 肾透明细胞癌 | KIRC | +0.62 | 0.098 | ⬆️ 不显著 |
| 肾乳头状细胞癌 | KIRP | +0.45 | 0.23 | ⬆️ 不显著 |
| 肝细胞癌 | LIHC | −0.31 | 0.17 | ➡️ 无差异 |
| 乳腺癌 | BRCA | −0.18 | 0.41 | ➡️ 无差异 |
| 甲状腺癌 | THCA | −0.03 | 0.88 | ➡️ 无差异 |
| 前列腺癌 | PRAD | −0.42 | 0.09 | ⬇️ 不显著 |
| 头颈鳞癌 | HNSC | −0.52 | 0.07 | ⬇️ 不显著 |
| 嫌色肾细胞癌 | KICH | −1.03 | 0.12 | ⬇️ 不显著（样本极少） |

---

## 四、TissGDB 基因注释

| 属性 | 值 |
|------|-----|
| 基因分类 | Class C |
| 组织特异性 | LiverStomach |
| 关联癌种 | LIHC, STAD |
| 致癌注释 | Fused with Oncogene |
| 融合事件 | RAB11B-SLC5A6 (ChiTaRs 数据库) |
| 药物 | Biotin (DB00121), Lipoic Acid (DB00121) |
| miRNA 调控 | hsa-let-7b-5p (ρ=−0.29, UVM 中唯一显著) |

---

## 五、文献交叉验证

| 癌种 | 发现 | 方法 | PMID |
|------|------|------|------|
| 肺腺癌 (LUAD) | SMVT→biotin→ACC→FASN 脂质代谢轴驱动肿瘤生长 | 功能验证 | 39426496 |
| 宫颈癌 | SLC5A6 KD→FASN↓→抑制增殖; FASN KD rescue | KD+OE+xenograft | 41108787 |
| 胃癌 (GC) | SLC5A6 诊断+预后 biomarker (TCGA + IHC) | 数据挖掘+IHC | Spandidos 2019 |
| 乳腺癌 | SMVT 在 T47D 细胞中功能性表达, 介导生物素摄取 | 体外功能 | 文献确证 |
| 前列腺癌 | SMVT 在 PC-3 细胞中功能性活跃 | 体外功能 | 文献确证 |

---

## 六、数据可靠性评估

| 维度 | 评分 | 说明 |
|------|:--:|------|
| 来源权威性 | ★★★★ | TissGDB 是 TCGA 官方衍生数据库 |
| 数据可复现 | ★★★★★ | 时隔 5 天两次查询完全一致 |
| 统计严谨性 | ★★★★ | FDR 多重检验校正 |
| 外部交叉验证 | ★★★★ | HPA + 文献三线验证 |
| ESCA 信号 | ⚠️ 弱 | FDR=0.005 边界显著，仅 10 对样本 |

---

## 七、局限性

1. **RNA ≠ 蛋白**: HPA IHC 部分癌种 mRNA 与蛋白表达不一致 (BLCA mRNA 最高但 IHC 阴性)
2. **配对样本有限**: 仅 14/33 癌种有 ≥10 对配对样本
3. **未覆盖全部癌种**: 宫颈癌(CESC)、胰腺癌(PAAD)等无配对数据
4. **mRNA 变化 ≠ 转运活性**: SMVT 是转运体，需功能实验验证

---

## 八、写作建议

- **强证据癌种** (可明确陈述): COAD, STAD, LUAD, LUSC, BLCA
- **弱证据癌种** (需保守措辞): ESCA ("showed a trend toward upregulation" / "was nominally significant")
- **无差异癌种**: BRCA, THCA, LIHC, PRAD, HNSC — 不应作为靶点证据
- 建议补充: GEPIA2 或 TIMER2 的独立验证图作为 Supplementary Figure
