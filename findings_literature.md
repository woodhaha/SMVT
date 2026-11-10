# SMVT 对接扩展新发现 — 文献意义调研

> 调研日期：2026-11-02 · 覆盖 2026-11-02 扩展筛选（11823 化合物）三大化学型新发现的文献定位

---

## 摘要（一句话）

这次扩展筛选最有价值的不是「多找了几个 hit」，而是**用计算手段把 SMVT 从「维生素转运体」重新框定为「抗惊厥药/中枢药的脱靶位点」**——既给 35 年已知的临床副作用（抗惊厥药致生物素缺乏）补上了分子机制，又提出了一个可合成的金刚烷羧酸先导骨架。

---

## 三大发现的意义（结合已发表文献）

### 1. 巴比妥/戊二酰亚胺 → 印证 1989 年临床发现（最强，最该写进稿件）

**关键文献**：Said HM, Redha R, Nylander W. *Biotin transport in the human intestine: inhibition by anticonvulsant drugs.* **Am J Clin Nutr 1989;49(1):127–31.** PMID 2911998 · DOI 10.1093/ajcn/49.1.127。

- 用肠刷状缘膜囊泡（BBMV）证明 **carbamazepine 和 primidone 竞争性抑制生物素转运**（Ki = 4.70 / 2.25 mmol/L），且对 D-glucose 转运无影响（特异性）。
- **primidone 在体内代谢为苯巴比妥** → 直接命中稿件 phenobarbital（ΔG −8.17）同一化学型。
- 后续临床研究：**>80% 长期服抗惊厥药的癫痫患者血浆生物素下降 45–50%**；phenytoin / phenobarbital / carbamazepine 均相关；抗惊厥药还加速生物素→bisnorbiotin 分解。

**意义**：稿件的 barbiturate 发现不是冷门计算巧合，而是**给一个 35 年已知的临床副作用（抗惊厥药致生物素缺乏）提供了分子层面的结构机制**。这是稿件「barbiturate 是 SMVT 非预期配体」最有力的外部背书。

**建议措辞**（Discussion）：
> "This docking finding recapitulates the long-standing clinical observation that anticonvulsants—including primidone, which is metabolised to phenobarbital—competitively inhibit SMVT-mediated intestinal biotin uptake and lower plasma biotin (Said et al., 1989; ...)."

### 2. 金刚烷/钻石烷羧酸 → 已知药化骨架，SMVT 连接是原创点

- 金刚烷羧酸本身是**成熟药化骨架**：`Discovery and optimization of adamantane carboxylic acid derivatives exhibiting potent diacylglycerol acyltransferase 1 (DGAT1) inhibitory activity`（PMID 26218650）。
- 即「金刚烷笼 + 羧基」chemotype 是 **drug-like、可合成、可优化**，不是凭空捏造。
- **但无任何文献把金刚烷羧酸与 SMVT 联系起来** → 这是稿件的**原创点**。

**意义**：稿件提出的金刚烷羧酸先导骨架既有药化可行性（DGAT1 抑制剂同骨架），又是全新的 SMVT 配体连接，可直接作为「新化学型」卖点。

### 3. 羧基药效团 → 与你自己的解释一致，已被实验坐实

- 已确立结论：**「SMVT 不与修饰过羧基的生物素类似物相互作用」**（biotin conjugates 综述，PMC10653662；另 PMID 16033951：desthiobiotin/泛酸强抑制生物素摄取，2-iminobiotin/diaminobiotin 弱）。
- 稿件的金刚烷羧酸簇（全带 COOH）符合此药效团；巴比妥无羧基，稿件已正确 hedge（"carboxyl-free scaffold, tentative"，因无 Na⁺ 对接会高估）。
- **本次检索证实这个 hedge 是对的**——羧基必需，巴比妥打分应作「暂定」。

### 4. NMDA 轴（MK-801/美金刚）→ 全新假设，但最弱

- **无文献**连接 NMDA 拮抗剂与 SMVT/生物素 → 纯原创假设。
- 但 MK-801 / memantine 是金刚烷**胺**（非羧酸），不符合「羧基锚定」药效团 → 对接分**最不可靠**，应作为最弱一档呈现，不宜过度强调。

---

## 关键参考文献清单

| # | 文献 | 要点 | 标识 |
|---|---|---|---|
| 1 | Said HM et al. Biotin transport in human intestine: inhibition by anticonvulsant drugs. *Am J Clin Nutr* 1989;49:127–31 | carbamazepine/primidone 竞争性抑制生物素转运 | PMID 2911998 |
| 2 | Biotin conjugates in targeted drug delivery (review) | SMVT 不与修饰羧基的生物素类似物作用 | PMC10653662 |
| 3 | Discovery/optimization of adamantane COOH DGAT1 inhibitors | 金刚烷羧酸 = 已知药化骨架 | PMID 26218650 |
| 4 | Structural basis for multivitamin recognition & transport by human SMVT | SMVT cryo-EM 结构（稿件已引 zhen2026） | s41467-026-74948-3 (Nat Commun 2026) |
| 5 | Foxd3/SLC5A6 axis regulates apoptosis in LUAD cells | SMVT 在肺腺癌的独立证据 | S0898656824004480 (2024) |
| 6 | Disturbances in biotin metabolism in children on long-term anticonvulsants | 抗惊厥药致生物素缺乏（临床） | Wiley 1998 |
| 7 | Sodium-dependent multivitamin transporter (SMVT): a potential target for drug delivery | SMVT 作为药物递送靶点（已知） | ResearchGate review |

---

## 对稿件的建议

> **✅ 已落地（2026-11-02）**：建议 1、2、4 已写入 `SMVT_manuscript.tex` Discussion + `references.bib`，编译 0 错误。建议 3（NMDA 降权）维持现状，未加文献。

1. **Discussion 补一段**：把「Said 1989 + 抗惊厥药生物素缺乏」作为 barbiturate 发现的临床印证 → ✅ 已加 `said1989`（line 428：carbamazepine/primidone 竞争性抑制 SMVT 生物素转运，Ki 4.70/2.25 mM）。
2. **金刚烷羧酸**：点一句 DGAT1 抑制剂同骨架 → ✅ 已加 `pagire2015`（line 440 结论）。
3. **NMDA 轴**：保持现状的轻描淡写（tentative），不加文献（无文献可加）。
4. **羧基 hedge**：补 cite 生物素类似物研究 → ✅ 已加 `park2005`（line 430，desthiobiotin/泛酸竞争 vs 2-iminobiotin/diaminobiotin 不竞争）。注：`tripathi2025`（已在稿）即 PMC10653662 那篇 biotin conjugates 羧基论文，故改用互补的 Park & Sinko 2005。

---

## 未解决 / 待验证

- [ ] 抗惊厥药抑制 SMVT 的**确切分子位点**（底物口袋竞争 vs 变构）——需 cryo-EM/突变验证，非本次 docking 可答
- [ ] 金刚烷羧酸是否**真被转运**（抑制 vs 底物）——需体外摄取实验（与 5-HTP 竞争摄取实验同设计）
- [ ] MK-801/memantine 的 SMVT 结合是否真实（无羧基，最不可靠）——建议降权或删
