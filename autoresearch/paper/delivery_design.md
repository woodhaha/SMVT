# COOH-bisCF3 偶联抗肿瘤药物递送设计 — 结直肠癌（CRC）及肝转移

> 2026-11-12 · 结合网络调研 + 计算研究（docking/RDKit）· 目标组织：CRC + CRC 肝转移

## 0. 背景与硬约束

- **SMVT–CRC 关联**：Neophytou 2022 (Cell Reports) — SMVT 在 Ras 驱动肠道肿瘤中诱导表达、沉默则缩瘤；SMVT 是肠道唯一 biotin 摄取系统，支持肠道干细胞有丝分裂。稿件自身数据：SMVT 在 COAD/READ 上调。
- **硬约束（已证）**：SMVT 介导的「分子偶联物转运」是死路——偶联物 22 Å 进不了 occluded 口袋（三态全正：AF2 +59.0 / 26va +23.7 / 26vd +3.2）。biotin 偶联物更小（618 Da）却因羧基被酯化耗尽、同样被拒（+14.2）。
- **可行机制只有两条**：纳米粒靶向胞吞（~100 nm 绕开口袋尺寸限制），或环糊精主客体。

---

## 方案 A — 纳米粒靶向胞吞（主推，最成熟）

```
COOH-bisCF3（SMVT 靶向头）─ 纳米粒（脂质体/PLGA）─ 装载 SN-38/伊立替康 或 奥沙利铂
        ↓ SMVT 识别（CRC 高表达）→ 介导胞吞 → 胞内释药
```

**靶向头**：COOH-bisCF3 = 最强含羧基 SMVT 配体（AF2 −9.68 / 26va −8.55 / 26vd −8.37，三态全第一），比 biotin（26va −6.49）强 ~2 kcal/mol。

**偶联化学（关键）**：用 **7 位（第 4 桥头）接 linker**，不动 1-羧基（SMVT 识别锚）和 3,5-CF3（亲和）。RQ5 已测 12 种 linker：7 位容忍**极小**取代基（甲基 −8.04），越大越糟（PEG/长链崩）。因此 linker 要短（单亚甲基），或直接走**环糊精⊃金刚烷主客体**（见方案 B 的几何限制）。

**Payload**（CRC 肝转移一线药）：
- **SN-38 / 伊立替康**（FOLFIRI 核心，伊立替康活性代谢物）——CRC 标准药。
- **奥沙利铂**（FOLFOX 核心，Pt 类）——有 adamantane-Pt(IV) 先例。

**先例**：biotin-SMVT 脂质体/cubosome 靶向（游离 biotin 竞争证实）；伊立替康纳米粒靶向 CRC（Baião 2026）；transporter-targeted nano-vehicles（Kou 2020）。

**优势**：机制最扎实（SMVT 胞吞有 biotin 先例，非转运）、临床相关性强、无 CD 腔几何约束。

---

## 方案 B — 分子级主客体（有 Pt 先例，但 CF3 撞 CD）

```
COOH-bisCF3（靶向头）─ 金刚烷⊃β-环糊精 ─ Pt(IV) 前药（奥沙利铂/顺铂）
```

**先例**：Shi 2012 (J Inorg Biochem) / Prashar 2011 (BMCL) — **adamantane-Pt(IV) 顺铂 + β-CD 主客体**；Ghitman 2023 adamantyl 偶联前药；Lee 2018 β-CD-PEG-叶酸-金刚烷靶向成像。

**⚠️ 几何冲突（本次 RDKit 发现）**：COOH-bisCF3 的 3,5-双 CF3 使分子最大核距 **7.3 Å**（金刚烷 3.5 Å，β-CD 内腔 ~6.0–6.5 Å）——即 vdW 尺度 ~10 Å，**远超 β-CD 腔**，CF3 会挡 β-CD 包合金刚烷笼。这使「COOH-bisCF3 ⊃ β-CD」在几何上不可行。

**三种修正（若坚持主客体）**：
1. **无 CF3 的普通金刚烷羧酸**（−7.73）做 guest——牺牲 CF3 亲和增益，恢复 β-CD 包合。
2. **双金刚烷设计**：COOH-bisCF3（靶向 SMVT）— linker — 普通金刚烷（⊃β-CD 装药）——两个金刚烷各司其职，但分子变大。
3. **共价 linker 替代主客体**——回到方案 A 的纳米粒/前药路线，避免 CD 腔约束。

---

## Payload 选择依据

| 药 | 机制 | CRC 地位 | 偶联可行性 |
|---|---|---|---|
| SN-38（伊立替康） | 拓扑异构酶 I 抑制剂 | FOLFIRI 核心，mCRC 一线 | C10 酚羟基酯偶联（伊立替康前药先例） |
| 奥沙利铂 | Pt-DNA 交联 | FOLFOX 核心，mCRC 一线 | Pt(II)→Pt(IV) 前药氧化偶联（adamantane-Pt 有先例） |
| 5-FU | 抗代谢 | FOLFOX/FOLFIRI 基础 | 太小，单用偶联意义弱 |

---

## 定稿：走递送胞吞路线（方案 A 为最终路线）

**决定**：采用**纳米粒靶向胞吞**，放弃 SMVT 转运 + β-CD 主客体（CF3-CD 几何冲突）。

### 具体设计产物（RDKit 已验证）

**靶向头（3 个偶联柄砌块）**：

| 名称 | 分子式 | MW | SMILES |
|---|---|---|---|
| 7-氨甲基-3,5-双(CF3)金刚烷-1-羧酸 | C₁₄H₁₇F₆NO₂ | 345.3 | `OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CN)` |
| 7-炔丙基-3,5-双(CF3)金刚烷-1-羧酸 | C₁₆H₁₆F₆O₂ | 354.3 | `OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CC#C)` |
| 7-叠氮甲基-3,5-双(CF3)金刚烷-1-羧酸 | C₁₄H₁₅F₆N₃O₂ | 371.3 | `OC(=O)C12CC3(C(F)(F)F)CC(CC(C3)(C1)C(F)(F)F)C2(CN=[N+]=[N-])` |

**完整偶联物（纳米粒 payload）**：

| 名称 | 分子式 | MW | 键型 |
|---|---|---|---|
| COOH-bisCF3-SN38 氨基甲酸酯前药 | C₃₇H₃₅F₆N₃O₈ | 763.7 | SN-38 C10 酚 → 氨基甲酸酯（pH/酯酶裂解）|
| COOH-bisCF3-奥沙利铂 Pt(IV) 前药 | 无机配合物 | — | Pt(II)→Pt(IV) 轴向羧酸酯 |

### 最终递送系统

```
COOH-bisCF3-NH2（靶向头, MW 345）
   └─ 酰胺键 ─ 脂质体/PLGA 纳米粒（~100 nm，表面多价展示靶向头）
        └─ 装载：SN-38（C10 酯）或 奥沙利铂 Pt(IV) 前药
   ↓ SMVT 识别（CRC 高表达）→ 介导胞吞 → 胞内/溶酶体释药
```

### 下一步验证
- 靶向头合成（市售 3,5-双(三氟甲基)金刚烷-1-羧酸为原料，7 位氨甲基化 1-2 步）。
- 纳米粒组装 + SMVT 高/低表达 CRC 细胞摄取对比 + 游离 biotin 竞争。

## 方案 A 合成路线（偶联化学）

### SN-38（伊立替康活性代谢物）
- SN-38 比伊立替康毒 200–2000×；伊立替康是 SN-38 的前药（羧酸酯酶水解）。
- **C10 酚羟基**是标准偶联位点 → **酯键**或**碳酸酯键**接纳米粒/前药（如 SN38-PA = 棕榈酸 C10 酯，亲脂前药自组装纳米粒）。
- 先例：Du 2017（双 7-乙基-10-羟基喜树碱偶联）、Qi 2024（IJN SN38 递送综述）、Wang 2025（SN38 前药纳米组装 CRC）。

### 奥沙利铂 Pt(IV) 前药
- Pt(II) 奥沙利铂 → **H₂O₂ 氧化** → Pt(IV) 增加两个**轴向配体**。
- 轴向羟基偶联：**氨基 → 异氰酸酯/氨基甲酸酯**；**羟基 → 碳酸酯（disuccinimidyl carbonate）**。
- 轴向配体可做靶向（葡萄糖 C50、alendronate、肽 R8K），也可做 linker 接纳米粒/抗体。
- 先例：Marotta 2023（Pt(IV) 前药综述）、O'Dowd 2023（奥沙利铂概述）、Sadeghi 2025。

### 靶向头偶联
COOH-bisCF3 通过 **7 位短 linker**（甲基/单亚甲基，RQ5 已测 −8.04）共价接纳米粒表面（酯/酰胺/点击），保留 1-羧基（SMVT 锚）+ 3,5-CF3（亲和）。

## 参考文献（已逐条 Crossref 核实，13/14 确认，DOI 已核）

1. Neophytou C, Pitsouli C. Biotin controls intestinal stem cell mitosis and host-microbiome interactions. *Cell Rep.* 2022. **10.1016/j.celrep.2022.110505** ✓
2. Shi Y, Dabrowiak JC. Host–guest interactions involving platinum anticancer agents. *Inorganica Chimica Acta.* 2012. **10.1016/j.ica.2012.06.021** ✓（期刊修正：非 J Inorg Biochem）
3. Prashar D, et al. Adamantane–platinum conjugate hosted in β-cyclodextrin. *Bioorg Med Chem Lett.* 2011. **10.1016/j.bmcl.2011.10.006** ✓
4. Hong Y, et al. Current status of nanoscale drug delivery systems for colorectal cancer liver metastasis. *Biomed Pharmacother.* 2019. **10.1016/j.biopha.2019.108764** ✓
5. Baião A, Castro, Viegas, et al. Targeted nanoparticle delivery of irinotecan. *J Control Release.* 2026. **10.1016/j.jconrel.2026.114609** ✓
6. Kou L, et al. Transporter-targeted nano-sized vehicles. *Cancers.* 2020. **10.3390/cancers12102837** ✓
7. Patel M, et al. SMVT molecular expression and functional activity (prostate). *Int J Pharm.* 2012. **10.1016/j.ijpharm.2012.06.011** ✓
8. Ghitman J, Voicu SI. Controlled drug delivery mediated by cyclodextrin-based supramolecular self-assemblies. *Carbohydr Polym Technol Appl.* 2023. **10.1016/j.carpta.2022.100266** ✓
9. Lee DW, et al. Supramolecular β-CD–adamantane for targeted cancer imaging. *J Ind Eng Chem.* 2018. **10.1016/j.jiec.2017.08.005** ✓
10. De Falco V, et al. How we treat metastatic colorectal cancer. *ESMO Open.* 2020. **10.1136/esmoopen-2020-000813** ✓（年份修正：2020 非 2019）
11. ~~Du Y 2017~~ — ⚠️ 标题未精确匹配，SN38 双偶联引用可用已核实的 Qi 2024（#12）替代。
12. Qi Q, et al. SN38 drug delivery system in cancer treatment. *Int J Nanomedicine.* 2024. **10.2147/IJN.S435407** ✓
13. Marotta C, et al. Recent advancements in anticancer Pt(IV) prodrugs. *Inorganica Chimica Acta.* 2023. **10.1016/j.ica.2023.121388** ✓
14. O'Dowd PD, et al. Oxaliplatin and its derivatives — an overview. *Coord Chem Rev.* 2023. **10.1016/j.ccr.2023.215439** ✓

> 核验方法：Crossref `query.bibliographic` 标题检索 + 作者/年份/期刊字段比对。13/14 确认真实，1 条（Du 2017）标题模糊未精确匹配，建议用已核实的 Qi 2024 替代或补精确标题后再核。
