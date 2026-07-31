# 5-HTP–SMVT 假说实验验证方案

> 创建: 2026-07-31 · 目的: 验证稿件 SMVT–5-HTP–serotonin 假说
> 科学问题: SMVT 是否介导 5-HTP 的 BBB 转运？NSAID 抑制是否影响脑 5-HT 合成？
> 依据: 稿件 docking (5-HTP pose 与 biotin 重叠 0.37Å) + Uchida 2015 (6 NSAID 抑制 SMVT)

---

## 计算预验证（2026-07-31）

5-HTP 具备 SMVT 底物识别的理化特征（支持 H1 前提）：

| 分子 | MW | logP | HBD | HBA | TPSA | COOH |
|------|----|----|-----|-----|------|------|
| Biotin | 244 | 0.80 | 3 | 3 | 78 | 1 |
| Pantothenic Acid | 219 | -1.04 | 4 | 4 | 107 | 3 |
| Lipoic Acid | 206 | 2.79 | 1 | 3 | 37 | 1 |
| **5-HTP** | **220** | **0.83** | **4** | **3** | **99** | **2** |

- **羧酸**: 5-HTP 有 α-氨基羧酸 → SMVT Na⁺-耦合羧酸识别位点
- **MW/logP 接近生物素**: 两亲性匹配底物识别
- **Tanimoto**: 5-HTP vs SMVT 底物 0.09-0.11（但 SMVT 底物彼此也只有 0.09-0.11）→ SMVT 识别特征而非结构，5-HTP 特征组合符合
- **Tryptophan vs 5-HTP = 0.675** → 结构相似，两者都可能走 LAT1，实验需区分 SMVT vs LAT1 贡献

→ 计算支持 5-HTP 为 SMVT 底物候选，需 Phase 1 摄取实验功能验证。

---

## 假说分解为 3 个可测命题

| # | 命题 | 当前证据 | 需验证 |
|---|------|---------|--------|
| H1 | 5-HTP 是 SMVT 底物 | docking 重叠 0.37Å（计算） | 摄取实验（功能） |
| H2 | NSAID 抑制 SMVT 介导的 5-HTP 摄取 | Uchida 2015（生物素/泛酸） | 5-HTP 特异性抑制 |
| H3 | 该抑制在 BBB 水平影响脑 5-HT | 无直接证据 | 体外 BBB + 体内 |

## 实验设计（按可行性排序）

### Phase 1: 体外摄取（最关键，验证 H1/H2）
**细胞系**: Caco-2（人肠上皮，内源 SMVT 高表达）或 SMVT-过表达 HEK293
**方法**: [3H]biotin 摄取作为阳性对照（已有 protocol），换 [3H]5-HTP 或冷 5-HTP + LC-MS/MS
**关键实验**:
1. **底物验证**: [3H]5-HTP 摄取是否 Na⁺ 依赖（换 Li⁺ 缓冲液）、饱和（Km 测定）、可被过量生物素/泛酸竞争抑制
2. **SMVT 特异性**: SMVT siRNA 敲低后 5-HTP 摄取下降幅度（对照: Uchida 2015 中 biotin 88.7%）
3. **NSAID 抑制**: diclofenac/indomethacin/ibuprofen 对 5-HTP 摄取的 IC50
**判定标准**: 5-HTP 摄取 Na⁺ 依赖 + 可被 SMVT 底物竞争 + siRNA 敏感 → H1 成立

### Phase 2: 体外 BBB 模型（验证 H3 的转运端）
**模型**: hCMEC/D3（Uchida 2015 用过的同一细胞系，已有 SLC5A6 表达数据）
**方法**: Transwell 双向转运实验，5-HTP 从管腔侧→脑侧通量
**关键实验**:
1. 5-HTP 跨膜转运是否被 SMVT 抑制剂（生物素过量 / NSAID）降低
2. 与 SLC5A6 表达量（WB/qPCR）的相关性

### Phase 3: 体内（验证 H3 整体）
**模型**: 小鼠，慢性 NSAID 给药（diclofenac，临床相关剂量）
**读出**: 脑组织 5-HTP/5-HT/5-HIAA 水平（HPLC）+ 行为学（可选）
**设计**: ± SMVT 抑制剂（生物素过量或药理剂量 NSAID），检测脑 5-HT 是否下降

## 预期结果与假说判定

| 结果 | 结论 |
|------|------|
| 5-HTP 摄取 Na⁺ 依赖 + siRNA 敏感 | SMVT 是 5-HTP 转运体（新发现） |
| NSAID 抑制 5-HTP 摄取 (IC50 合理) | NSAID-SMVT-5-HTP 抑制轴成立 |
| 体内 NSAID 降脑 5-HT | 完整机制链条 |

## 关键对照
- **阴性对照**: 5-HTP 摄取的 Na⁺ 非依赖部分（可能被动扩散）、非转运蛋白对照
- **阳性对照**: biotin 摄取（已知 SMVT 底物）、LAT1 抑制剂 (BCH) 区分 5-HTP 是否也走 LAT1
- **特异性**: 5-HTP vs 色氨酸（可能都走 LAT1，需区分 SMVT 贡献）

## 与稿件的关系
- 若 H1 证实：稿件假说从"推测"升级为"计算 + 实验支持"，投稿后 rebuttal 有底牌
- 稿件 Limitations 已声明 "5-HTP competition uptake assays in SMVT-expressing cell lines and BBB models are needed" — 此方案即对应

## 资源需求
- [3H]5-HTP（PerkinElmer 定制）或冷 5-HTP + LC-MS
- hCMEC/D3（已有文献protocol，Uchida 实验室方法可参考）
- 小鼠 + HPLC（温医一院药学平台可用性待查）
