# SLC5A6 (SMVT) — TCGA 泛癌差异表达分析

> 查询日期: 2026-06-23
> 基因: SLC5A6 / SMVT (Sodium-dependent Multivitamin Transporter)
> Entrez ID: 8884 | Cytoband: 2p23.1

---

## 1. Tumor vs Normal 显著差异癌种

来源: TissGDB (TCGA IlluminaHiSeq_RNASeqV2, pan-cancer normalized log2(norm_counts+1))

| 癌种 | 缩写 | Tumor均值 | Normal均值 | Log2FC | P-value | FDR | 趋势 |
|------|------|-----------|------------|--------|---------|-----|------|
| 膀胱尿路上皮癌 | BLCA | 2.23 | 0.52 | **+1.71** | 2.40×10⁻⁸ | 1.84×10⁻⁶ | ⬆️ Up |
| 肺鳞状细胞癌 | LUSC | 1.96 | 0.28 | **+1.68** | 1.29×10⁻²⁰ | 1.64×10⁻¹⁹ | ⬆️ Up |
| 结肠腺癌 | COAD | 2.78 | 1.27 | **+1.51** | 2.09×10⁻¹² | 4.87×10⁻¹¹ | ⬆️ Up |
| 食管癌 | ESCA | 2.11 | 0.63 | **+1.48** | 1.87×10⁻⁴ | 5.01×10⁻³ | ⬆️ Up |
| 胃腺癌 | STAD | 1.96 | 0.52 | **+1.44** | 4.68×10⁻⁹ | 2.77×10⁻⁷ | ⬆️ Up |
| 肺腺癌 | LUAD | 1.29 | 0.16 | **+1.13** | 8.99×10⁻²⁶ | 4.20×10⁻²⁴ | ⬆️ Up |

> 仅有 ≥10 对 Tumor-Normal 配对样本的癌种被纳入差异分析 (共14种)。
> Log2FC = log2(Tumor_mean) - log2(Normal_mean)

**关键结论**: SLC5A6 在这6种癌种中**一致上调**，尤以肺鳞癌(Log2FC=1.68, P=1.3×10⁻²⁰)、结肠癌(Log2FC=1.51, P=2.1×10⁻¹²)最为显著。

---

## 2. 泛癌 RNA 表达分类 (Human Protein Atlas)

| 属性 | 值 |
|------|-----|
| RNA 组织特异性 | **Expressed in all** — 在全部17种TCGA癌种中均有表达 |
| RNA 癌症特异性分类 | Low cancer specificity (泛表达) |
| 蛋白表达定位 | 细胞质 + 膜 (Cytoplasmic and membranous) |
| IHC 蛋白表达 | 鳞状细胞癌、结直肠癌、胰腺癌、尿路上皮癌呈弱-中等阳性，其余肿瘤阴性 |
| 正常组织蛋白表达 | 多种组织中细胞质和膜表达 |

---

## 3. 已知促癌机制

```
SLC5A6 → 生物素摄取↑ → ACC ↑ → FASN ↑ → 脂质合成↑ → 肿瘤增殖
              │
              └──→ 辅酶A合成 → 线粒体能量代谢
```

### 文献证据

| 癌种 | 发现 | 实验方法 | 证据等级 | 文献 |
|------|------|---------|---------|------|
| 肺腺癌 (LUAD) | SMVT→biotin→ACC→FASN脂质代谢轴驱动肿瘤生长 | 功能验证 | ★★★ | PMID: 39426496 |
| 宫颈癌 | SLC5A6敲低→FASN↓→抑制增殖; FASN KD可逆转SLC5A6 OE效应 | KD+OE+rescue+xenograft | ★★★ | PMID: 41108787 |
| 胃癌 (GC) | SLC5A6 诊断+预后biomarker (TCGA验证) | TCGA数据挖掘+IHC | ★★☆ | Spandidos, 2019 |
| 乳腺癌 | SMVT在T47D细胞中功能性表达,介导生物素摄取 | 体外功能验证 | ★★☆ | 文献确证 |
| 前列腺癌 | SMVT在PC-3细胞中功能性活跃,介导生物素摄取 | 体外功能验证 | ★★☆ | 文献确证 |

---

## 4. 对 SMVT 靶点研究的启示

### 4.1 作为药物递送靶点
- **优势**: SMVT在肿瘤中高表达,特别是在消化道(COAD/STAD/ESCA)和肺部(LUSC/LUAD)
- **策略**: 利用SMVT底物(生物素/泛酸/硫辛酸)进行**靶向药物偶联**(SMVT-targeted drug conjugate)
- **已知底物**: Biotin (Km~), Pantothenate, Lipoic acid, Iodide

### 4.2 作为治疗靶点
- **机制**: 抑制SMVT→切断肿瘤细胞的生物素供应→阻断ACC-FASN脂质合成通路
- **注意事项**: 正常组织(肠道、肾脏)也表达SMVT → 需评估治疗窗口
- **选择性窗口**: 肿瘤细胞对脂质合成的依赖性远高于正常细胞("metabolic addiction")

### 4.3 诊断/预后价值
- 胃癌: TCGA mRNA + IHC蛋白双验证,潜在诊断biomarker
- 6种癌种中Tumor vs Normal均有显著差异 → 泛癌诊断潜力

### 4.4 关键未解决问题
1. SLC5A6 上调的**转录调控机制**是什么?(哪个转录因子驱动?)
2. 为什么某些癌种(如BLCA)表达最高但IHC反而阴性?
3. SMVT介导的药物偶联物**肿瘤选择性**如何?
4. SMVT在**免疫微环境**中的角色?

---

## 5. 数据来源

| 来源 | URL |
|------|-----|
| TissGDB | https://bioinfo.uth.edu/TissGDB/ (SLC5A6, Tissue ID: 8884) |
| Human Protein Atlas | https://www.proteinatlas.org/ENSG00000138074-SLC5A6 |
| NCBI Gene | https://www.ncbi.nlm.nih.gov/gene/8884 |
| Affinage | https://affinage.wi.mit.edu/gene/SLC5A6 |
| PMC (SLC family review) | https://pmc.ncbi.nlm.nih.gov/articles/PMC10308049/ |
