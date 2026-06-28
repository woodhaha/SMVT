# SMVT (SLC5A6) — 靶点研究报告

> 生成日期: 2026-06-23  
> 数据来源: DrugBank, PubMed, OpenAlex, ClinicalTrials.gov  
> 工具: DrugBank-Database Skill + WebSearch

---

## 1. 基础信息

| 项目 | 详情 |
|------|------|
| **基因** | *SLC5A6* |
| **蛋白** | Na⁺-依赖多维生素转运体 (SMVT) |
| **家族** | SLC5 Solute:Sodium Symporter Family |
| **UniProt** | Q9Y289 |
| **GeneCards** | [SLC5A6](https://www.genecards.org/card/SLC5A6) |
| **转运机制** | Na⁺ 偶联同向转运 (2 Na⁺ : 1 底物) |

---

## 2. 生理底物

| 底物 | 类型 | 亲和力 |
|------|------|--------|
| **Biotin (维生素 B₇)** | 主要底物 | 高亲和力 |
| **Pantothenic Acid (维生素 B₅)** | 主要底物 | 中等亲和力 |
| **α-Lipoic Acid (R-对映体选择性)** | 底物 | 中等亲和力 |
| **Iodide (I⁻)** | 意外底物 (2011 年发现) | 低亲和力 |
| **Na⁺** | 共转运离子（驱动力） | 必需 |

---

## 3. DrugBank 收录的药物与底物

| 药物/底物 | DrugBank ID | 状态 | 关系 |
|-----------|-------------|------|------|
| Biotin (维生素 B₇) | DB00121 | Approved / Nutraceutical | 内源性底物 |
| α-Lipoic Acid | DB00166 | Approved / Nutraceutical | 内源性底物 |
| Pantothenic Acid (B₅) | — | Nutraceutical | 内源性底物 |
| Iodide (I⁻) | — | 离子 | 底物 (de Carvalho & Quick, *J Biol Chem*, 2011) |
| **Gabapentin enacarbil (Horizant®)** | — | **FDA Approved** | **唯一获批 SMVT 靶向药物** — 转运前药 |

> Gabapentin enacarbil 利用肠道 SMVT 提高加巴喷丁口服生物利用度，适应症：带状疱疹后神经痛。

---

## 4. SMVT 药理学抑制

### 竞争性底物 (相互抑制)

Biotin、Pantothenic acid、Lipoic acid 竞争同一结合位点，相互抑制转运。Iodide 抑制所有有机底物的摄取。

### 外源性抑制剂 (BBB — hCMEC/D3 验证)

| 类别 | 化合物 |
|------|--------|
| **NSAIDs** | Indomethacin、Ketoprofen、Diclofenac、Ibuprofen、Phenylbutazone、Flurbiprofen |
| **脂质介质** | Prostaglandin E₂ (PGE₂) |
| **脂肪酸** | Docosahexaenoic acid (DHA) |

> 来源: Uchida et al. (2015), *Journal of Neurochemistry*

---

## 5. SMVT 激活剂

| 化合物 | 机制 |
|--------|------|
| D-(+)-Biotin | 直接底物结合 → 驱动转运 |
| α-Lipoic acid | 直接底物 → 增强功能活性 |
| D-Panthenol | 泛酸结构类似物 |
| Biotin-PEG-amine 衍生物 | 与 Biotin 结构相似 |
| **Amiloride** | 抑制 Na⁺ 通道 → 增加 Na⁺ 梯度驱动力 |
| **Ouabain / Digoxin** | 抑制 Na⁺/K⁺ ATPase → 升高胞内 Na⁺ |
| **Bumetanide** | 抑制 NKCC → 升高胞内 Na⁺ |

---

## 6. 组织分布与治疗意义

| 组织 | 表达 | 治疗意义 |
|------|------|----------|
| 小肠上皮 | 高 | **口服前药递送** (如 gabapentin enacarbil) |
| 肝细胞 | 高 | 肝脏靶向 |
| 肾近端小管 | 高 | 肾重吸收 / 肾毒性规避 |
| 血脑屏障 (微血管内皮) | 中-高 | CNS 药物靶向 |
| 乳腺 | 中 | 乳腺癌靶向 |
| **乳腺肿瘤** | **过表达** | **肿瘤选择性靶向** |
| **前列腺肿瘤** | **过表达** | **肿瘤选择性靶向** |
| **视网膜母细胞瘤** | **过表达** | **肿瘤选择性靶向** |
| **胃癌** | **过表达** | **肿瘤选择性靶向** |
| **胶质瘤** | **过表达** | **肿瘤选择性靶向** |
| 角膜/视网膜上皮 | 中 | 眼部给药 |
| 胎盘 | 有表达 | 母胎营养转运 |

---

## 7. SMVT 靶向药物递送管线 (2024–2025)

### 7.1 已获批

| 药物 | 适应症 | 机制 |
|------|--------|------|
| Gabapentin enacarbil (Horizant®) | 带状疱疹后神经痛 | SMVT 转运前药 |

### 7.2 临床阶段

| 疗法 | 临床试验 | 适应症 | 状态 |
|------|----------|--------|------|
| [⁶⁸Ga]Ga-FAPI-Biotin PET/CT | NCT06740240 | 多种癌症诊断 | 启动中 (2024.12) |
| [¹⁷⁷Lu]DOTA-Biotin / AvidinOX | — | 结直肠癌肝转移 | 剂量递增 |

### 7.3 临床前活跃管线

| 策略 | 载荷 | 文献 |
|------|------|------|
| **Biotin-小分子偶联物** | Ursolic acid | *Molecules*, 2025 |
| | Resorcinol-PD-L1 抑制剂 | *J Med Chem*, 2023 |
| | Podophyllotoxin 衍生物 | — |
| | 铂(IV) 复合物 | 乳腺癌 |
| | Paclitaxel (延长 spacer) | — |
| | Gemcitabine-Coumarin (theranostic) | — |
| **Biotin-纳米粒** | PAMAM G4.5 dendrimer | 多种 |
| | 固体脂质纳米粒 (MTX) | *ACS Appl Nano Mater*, 2025 |
| | 热响应性胶束 | — |
| | Quercetin + Doxorubicin 共包载 | 乳腺癌耐药逆转 |
| | 谷胱甘肽双靶向 (BBB 穿透) | CNS |
| | 胆红素纳米粒 (TME 响应) | — |

---

## 8. SMVT 调控

### 转录调控

- **KLF-4**、**AP-2** 转录因子
- 全羧化酶合成酶 (HCS) 介导的染色质重塑

### 翻译后调控

- **PKC**、**Ca²⁺/钙调蛋白**、**酪蛋白激酶 2** 通路
- N-糖基化
- PDZ 结构域蛋白互作 (**PDZD11**)

### 炎症下调

- LPS 和促炎细胞因子 → 降低 SMVT 膜表达

---

## 9. 临床遗传学

- **SLC5A6 功能丧失突变** → 钠依赖型多维生素转运体缺乏症
- 表型：发育迟缓、神经病变、免疫缺陷
- **可治疗**：大剂量 Biotin + Pantothenate + Lipoate 补充
- 2024 年新发现的中间表型变体扩展了表型谱 (*J Hum Genet*, 2024)

---

## 10. 2024–2025 关键争议与进展

| # | 进展 | 来源 |
|----|------|------|
| 1 | **Biotin-纳米粒是否全部经 SMVT 进入细胞？** MCTs 和内吞途径可能共同参与 | *IJMS*, Feb 2025 |
| 2 | *Slc5a6* 敲除小鼠揭示代谢性扩张型心肌病 → 维生素干预可挽救 | bioRxiv, June 2025 |
| 3 | *SLC5A6* 基因失活测试系统建立 | *Acta Naturae*, Oct 2025 |
| 4 | 综合综述：Biotin 抗癌治疗 | *Nanoscale*, Jan 2025 |
| 5 | SMVT 与受体/内吞途径在 Biotin 偶联物内化中的角色再评估 | *J Enzyme Inhib Med Chem*, 2023 & 2025 |

---

## 11. 研究工具与资源

| 资源 | 链接 |
|------|------|
| DrugBank | https://go.drugbank.com |
| VARIDT Transporter DB | https://varidt.idrblab.net/data/transporter/details/dtd0426 |
| GeneCards SLC5A6 | https://www.genecards.org/card/SLC5A6 |
| Santa Cruz SMVT Activators | https://www.scbio.cn/browse/smvt_slc5a6-activators |
| ClinicalTrials.gov ⁶⁸Ga-FAPI-Biotin | NCT06740240 |

---

## 12. 核心结论

1. **DrugBank 直接收录的 SMVT 底物有限**（3 个内源性 + 1 个前药），但 SMVT 的**肿瘤过表达**特性使其成为极其活跃的药物递送靶点
2. **Gabapentin enacarbil** 是唯一上市药物，验证了 SMVT 靶向递送的临床可行性
3. **Biotin 偶联策略**是 2024–2025 最热门的方向——小分子、纳米粒、放射治疗三位一体
4. **核心未解问题**：Biotin 偶联物的内化机制是否完全依赖 SMVT，还是 MCTs + 内吞共同参与？
5. **NSAIDs 对 SMVT 的抑制**值得关注——临床 NSAID 用药可能干扰 SMVT 靶向递送
