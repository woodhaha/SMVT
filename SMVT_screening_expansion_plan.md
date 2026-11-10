# SMVT 虚拟筛选库扩展方案

> 创建：2026-10-30 · 目的：把对接库从 ~703 化合物扩到万级，强化「SMVT 是非预期药物转运底物 → DDI 风险」叙事 + 找新化学型

---

## 1. 背景与目标

SMVT (SLC5A6) 是 Na⁺ 依赖的多维生素转运体（底物：生物素、泛酸、硫辛酸）。当前稿件已用 **703 化合物**（FDA + ChEMBL-FDA + DrugBank top500）对接出 **8 个 elite hits**（ΔG < −8 kcal/mol），发现巴比妥类/NSAID/阿片类为「非预期 SMVT 配体」。

**扩展目标**：
1. 广覆盖 —— 从「已批准药」扩到「临床期 + 在研 + 生物活性 + 天然产物」
2. 强叙事 —— 更多药物类别命中，坐实 DDI 风险
3. 找新骨架 —— 超越已知药，引入维生素/氨基酸类似物化学型

---

## 2. 现状盘点

| 来源 | 数量 | 文件 |
|---|---|---|
| FDA 手挑 | 54 | `03_Analysis/data/fda_drugs.csv` |
| ChEMBL FDA 批准 | 200 | `03_Analysis/data/chembl_fda_approved.csv` |
| DrugBank top500 | 500 | `03_Analysis/outputs/drugbank_top500_for_docking.csv` |
| **合计** | **~703** | `03_Analysis/outputs/screening_master_results_ALL.csv` |

对接配置：`03_Analysis/data/docking_config.json`（中心 5.03, 1.5, −2.44；盒子 25³；exhaustiveness 16）。

---

## 3. 候选库清单（按优先级）

| # | 库 | 规模 | 开放度 | 格式 | 对本稿价值 |
|---|---|---|---|---|---|
| ① | **Drug Repurposing Hub (Broad)** | **7,423** 化合物 | ✅ CC-BY 4.0 免费 | SMILES + 靶点/MOA/临床期 | 临床期 + 在研，直接撑 DDI 叙事 |
| ② | **ChEMBL 全生物活性**（SMVT 相关靶点 + drug-like 过滤） | ~2.4M（可筛到 ~10K–100K） | ✅ 完全开放 | SMILES/InChI | 广覆盖，已有 client 脚本 |
| ③ | **ZINC drug-like 子集** | 数百万（可采购） | ✅ 免费 | SMILES/SDF/mol2 | 新化学型，需抽样 |
| ④ | **COCONUT 天然产物** | ~40 万 | ✅ 免费 | SDF/CSV/SMI | 维生素/氨基酸类似物新骨架 |
| ⑤ | **DrugBank 全库** | ~11,000 | ⚠️ 学术账号 | XML/SDF | 你已有 2823 评分，可扩全 |
| ⑥ | **底物类似物枚举**（biotin/pantothenate/lipoate） | 数百~数千（生成） | 本地生成 | SMILES | **SMVT 特异性最强** |

---

## 4. 分阶段策略

### 阶段 A —— 立即可做（下载 + 对接，量级可控）
**Drug Repurposing Hub 全库 7,423 化合物**

- 下载（2025-08-18 版，已验证 URL）：
  - 注释（名称/临床期/MOA/靶点）：`https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20250818.txt`
  - 样品（含 **SMILES**/InChIKey/PubChem ID）：`https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20250818.txt`
- 许可：CC-BY 4.0（可引用，注明 Corsello 2017 Nat Med）
- 处理：合并两表 → 去盐/去重 → 过滤掉已对接的 703 → 得 ~6,800 新化合物
- 对接：Vina exhaustiveness 16，~6,800 × 1–2 min = **~115–225 CPU·小时**（8–16 并行实例可压到 1–2 天）
- **增量价值**：临床期/在研药是「已批准药」的直接延伸，命中可直接写 DDI 预警

### 阶段 B —— 广度扩展（需过滤/抽样）
**ChEMBL + ZINC**

- ChEMBL：用现有 `download_chembl_client.py`，查 SMVT/SLC5A6 相关靶点（SLC 家族、维生素转运、生物素/泛酸代谢）活性化合物 + 相似性搜索（以 8 个 hit + 3 个底物为 query）
- ZINC：下 `drug-like` 子集，按 Tanimoto 相似度对底物/hit 抽样 ~10K，避免盲目全筛
- 目标：新增 ~10K–20K 化合物

### 阶段 C —— 新骨架（针对 SMVT 底物）
**COCONUT 天然产物 + 底物类似物枚举**

- COCONUT：下 SDF（`coconut.naturalproducts.net/download` 或 Zenodo），过滤 MW<500 + 含羧基（SMVT 底物关键特征）+ 相似底物
- 底物类似物枚举：RDKit R-group 枚举围绕生物素（ureido 环/戊酸链/羧基）、泛酸（β-丙氨酸/泛解酸）、硫辛酸（二硫戊环/戊酸链）
- 目标：新增 ~1K–5K 高针对性化合物

---

## 5. 执行细节

每个库的处理流水线（复用现有 `docking_pipeline.py` 的 ligand 准备 + Vina）：

1. **下载** → 2. **RDKit 清洗**（去盐、标准化、质子化 pH 7.4）→ 3. **去重**（InChIKey，与现有 703 交叉去重）→ 4. **Lipinski/drug-like 过滤**（RO5 + 可旋转键 ≤10）→ 5. **3D 构象**（RDKit ETKDG + MMFF 优化）→ 6. **Vina 对接**（配置复用 `docking_config.json`）→ 7. **结果合并**（与现有 `screening_master_results_ALL.csv` 同 schema）

---

## 6. 计算成本估算

| 阶段 | 化合物数 | 单化合物耗时 | 总耗时 | 并行方案 |
|---|---|---|---|---|
| A（Repurposing Hub） | ~6,800 | 1–2 min | 115–225 CPU·h | 16 核并行 ≈ 8–14h |
| B（ChEMBL+ZINC） | ~10–20K | 1–2 min | 170–670 CPU·h | 32 核并行 ≈ 5–20h |
| C（COCONUT+枚举） | ~1–5K | 1–2 min | 17–170 CPU·h | 16 核 ≈ 1–10h |

⚠️ **资源冲突**：ComputeShare pod（RTX 4090）当前被 Na⁺ MD 占满（~25h）。对接是 CPU 密集，可用本地 CPU 或另开 pod；也可等 MD 跑完。

---

## 7. SMVT 特异性策略（最高命中率路径）

SMVT 底物识别关键特征（已从稿件 findings 提炼）：
- **羧基**（biotin/5-HTP 都有，巴比妥酸缺羧基 → 不能转运）
- MW ~220–250（biotin 244、5-HTP 220）
- 极性 H 供体/受体（氨基、羟基）

因此**最高价值不是盲目扩库，而是枚举含羧基的小分子类似物**：
- 生物素类似物：改造 ureido 环 / 戊酸侧链 / 羧基位置
- 泛酸类似物：β-丙氨酸 + 泛解酸组合
- 硫辛酸类似物：二硫戊环 + 戊酸
- 氨基酸类似物（5-HTP 同类）：色氨酸/苯丙氨酸衍生物

这条路径可直接用 RDKit 本地生成（零下载），且命中率远高于随机库。

---

## 8. 风险与对策

| 风险 | 对策 |
|---|---|
| 大库盲目对接命中率低 | 先做相似性/药效团预筛，再对接 top-N |
| 羧基缺失导致假阴性（SMVT 特异性） | 保留羧基/极性基过滤，不盲目套 RO5 |
| 计算资源冲突（MD 占 pod） | 对接走本地 CPU 或另开 pod，与 MD 并行 |
| 命中化合物非可采购 | 优先 Repurposing Hub（可采购）+ ZINC（标注可采购） |
| 引用/许可 | Repurposing Hub 用 CC-BY 4.0 需引 Corsello 2017；COCONUT 引 COCONUT 2.0 2025 |

---

## 9. 决策点（待你定）

1. **先跑阶段 A（Repurposing Hub 6,800）吗？** 还是 A+B+C 一起排期？
2. **对接资源**：本地 CPU 并行 / 另开 pod / 等 MD 跑完再用现 pod？
3. **底物类似物枚举**（§7）是否作为独立高优先路线先做？

---

## 附：关键 URL 速查

- Drug Repurposing Hub 注释：`https://repo-hub.broadinstitute.org/public/data/repo-drug-annotation-20250818.txt`
- Drug Repurposing Hub 样品（SMILES）：`https://repo-hub.broadinstitute.org/public/data/repo-sample-annotation-20250818.txt`
- ZINC 下载：`https://files.docking.org`（drug-like 子集）
- ChEMBL API：`https://www.ebi.ac.uk/chembl/api/data/`
- COCONUT 下载：`https://coconut.naturalproducts.net/download`
- NPASS：`http://bidd.group/NPASS/`
