# SMVT 论文代码归档计划（PLOS CB Data/Code Availability）

> PLOS CB 硬性要求：投稿时 code 必须公开在 repository（GitHub/Zenodo），不是 "upon acceptance"。

## Repository 结构（建议）

```
smvt-pancancer/
├── 01_expression_survival/     # 泛癌表达 + 生存
├── 02_mutation_esm2/           # 突变景观 + ESM-2 致病性
├── 03_ppi_pathway/             # PPI + 通路富集
├── 04_virtual_ko/              # 虚拟 KO
├── 05_virtual_screening/       # 虚拟筛选（ML + docking + SAR）
├── 06_md_mmgbsa/               # MD + MM-GBSA
├── README.md                   # 复现说明
└── data_sources.md             # 数据来源清单
```

## 核心脚本清单（结果依赖的，需公开）

### 01 表达 + 生存
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `run_pancancer.py` | DE 分析（Mann-Whitney + BH-FDR） | TCGA TOIL FPKM |
| `multivariate_survival.py` | Cox 多变量（age/sex/stage） | TCGA 临床 |
| `survival_analysis_real.py` | 生存分析（33 癌型） | TCGA 临床 |
| `survival_gdc.py` | GDC 生存数据提取 | NCI GDC API |
| `fig_pancan_de.py` | Fig1 泛癌表达 | 上列输出 |

### 02 突变 + ESM-2
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `extract_slc5a6_maf.py` | MC3 MAF 预提取 SLC5A6 | MC3 MAF（3.5G） |
| `fig_mutation_final.py` | Fig3 突变景观 | 提取的 CSV |
| `esm_pathogenicity_mc3.py` | ESM-2 35M 致病性 | HF esm2_t12_35M |
| `esm650_crossval.py` | ESM-2 650M 交叉验证 | HF esm2_t33_650M |

### 03 PPI + 通路
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `visualize_string_network.py` | STRING 网络（Fig5/6） | STRING API |
| `pathlinkR_analysis.R` | 通路富集（FigS4/S5） | pathlinkR R 包 |

### 04 虚拟 KO
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `run_scTenifoldKnk_light.R` | scTenifoldKnk KO（Table DRG） | GSE178341 scRNA |
| `pc_regression_GRN.py` | PC 回归验证（ρ=0.712） | 同上 |

### 05 虚拟筛选（核心，最多脚本）
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `pharmacophore_ml_screen.py` | ML 药效团筛选（AUC 0.888） | 84 训练化合物 |
| `docking_batch_screen.py` | 主 docking（500 化合物） | AF-Q9Y289-F1 |
| `docking_expanded.py` | 11,823 扩展库 docking | Repurposing Hub + Selleck |
| `merge_anticancer.py` | 扩展库合并去重 | 三源化合物库 |
| `pathc_v2.py` | 巴比妥/戊二酰亚胺聚焦 | PubChem 类似物 |
| `pathc_barbiturate_hunt.py` | 巴比妥 SAR | PubChem |
| `glutarimide_hunt.py` | 戊二酰亚胺 SAR | PubChem |
| `sar_deep_analysis.py` | 金刚烷 58 类似物 SAR | 本地 docking |

### 06 MD + MM-GBSA
| 脚本 | 用途 | 数据依赖 |
|------|------|---------|
| `md_v13.py`（⚠️ 云端 weste） | MD 协议（HMR + barostat） | OpenMM 8.5.2 |
| `SMVT-MD-Analysis/*.py` | MD 分析（RMSD/FEL/DCCM） | DCD 轨迹 |
| `MM-GBSA Analysis/*.py` | MM-GBSA 能量分解 | 轨迹 |

## ⚠️ 三个投稿前必须解决的缺口

1. **`md_v13.py` 在云端**（weste 服务器已关机）——需从云端同步回本地，或标注"MD 协议脚本备索"
2. **论文 Code availability 声明错误**——说 `MD_pipeline/`，实际 MD 脚本在 `SMVT-MD-Analysis/` + `MM-GBSA Analysis/`，需修正
3. **措辞错误**——当前 "deposited upon acceptance"，PLOS CB 要求投稿时已公开，需改成 "available at [GitHub/Zenodo URL]"

## 数据来源清单（data_sources.md 内容）

- TCGA（NCI GDC）：RNA-seq + 临床 + MC3 MAF
- GEO GSE178341：结直肠癌 scRNA-seq（371,000 细胞）
- gnomAD v4：pLI/LOEUF/%HI 约束指标
- STRING v12：PPI 网络
- ChEMBL v34：3,311 已批准药
- Drug Repurposing Hub（7,423）+ Selleck L3000 + WHO ATC L01 + NCI（扩展库）
- AlphaFold AF-Q9Y289-F1 + cryo-EM PDB 26va–26ve
- HF：esm2_t12_35M_UR50D + esm2_t33_650M_UR50D

## 执行顺序建议

1. 建 GitHub 私有 repo（投稿前转公开），按上述 6 模块归档脚本
2. 从 weste 同步 `md_v13.py`（或标注备索）
3. 写 README.md（复现步骤）+ data_sources.md
4. 修正论文 Code availability 段（路径 + "available at URL" 措辞）
5. Zenodo 存档（拿 DOI）或 GitHub release
