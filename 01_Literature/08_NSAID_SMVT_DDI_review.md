# NSAID–SMVT 药物相互作用综述（文献挖掘）

> 整理: 2026-07-31 · 目的: 支撑稿件 SMVT–5-HTP–serotonin 假说的机制链条
> 关联: 稿件 Discussion (5-HTP 段) + 虚拟筛选 NSAID 命中

---

## 机制链条

```
NSAID ──抑制──▶ SMVT (SLC5A6) ──介导──▶ BBB 生物素/泛酸/5-HTP 摄取
                │
                └──▶ 脑内 5-HTP 下降 → 5-HT 合成减少 → 抑郁风险↑
```

## 关键证据

### 1. SMVT 是 BBB 生物素/泛酸主要转运体（实验确证）
**Uchida 2015, *J Neurochem* 134:97–112, DOI 10.1111/jnc.13092, PMID 25809983**
- SLC5A6 表达于人类脑微血管，定位于脑毛细血管内皮**腔膜面**（luminal）
- hCMEC/D3 细胞中 SLC5A6 介导 **88.7% [3H]biotin / 98.6% [3H]pantothenic acid** 摄取（siRNA 敲低验证）
- **6 个 NSAID 显著抑制 SLC5A6 转运**：indomethacin, ketoprofen, diclofenac, ibuprofen, phenylbutazone, flurbiprofen
- 另被抑制物：PGE2, lipoic acid, DHA
- 意义：首次确认人类 BBB 中 SLC5A6 是生物素/泛酸 luminal 摄取的主要贡献者

### 2. 虚拟筛选复现 NSAID 结合（本工作计算验证）
- Diclofenac −7.15, flufenamic −6.71, niflumic −6.68, meclofenamic −6.63, tolfenamic −6.10 kcal/mol（Vina）
- 与 Uchida 实验抑制的 NSAID 化学型（乙酸类/丙酸类/芬那酸类）重叠
- 5-HTP docked pose 与 biotin 几乎重叠（质心差 0.37Å, FigS15）→ SMVT 底物位点可容纳 5-HTP

### 3. NSAID 使用与抑郁风险（流行病学，证据复杂）
| 研究 | 人群 | 结果 |
|------|------|------|
| Molero 2019 (SUN 队列) | 22,564 一般人群 | 常规 NSAID 使用与抑郁**无关联** |
| Hu 2020 (BMC Med, 瑞典注册) | 316,904 癌症患者 | **非阿司匹林 NSAID 增加抑郁风险 HR 1.24**；COX-2 抑制剂 HR 1.26；阿司匹林保护性 HR 0.88 |
| Shaikh 2023 (Cancer Inform, DOI 10.1177/11769351231165161) | 14,992 老年癌症+OA 幸存者 | 累积 NSAID 天数与抑郁正相关（第 6 大预测因子），90/120 天时抑郁率最高 25% |

**解读**：一般人群无显著关联，但癌症/老年高风险人群中非阿司匹林 NSAID 与抑郁风险正相关——与 SMVT 抑制机制一致（癌症患者 SMVT 表达高，抑制效应更显著）。

### 4. 5-HTP 的 BBB 转运空白（假说的新颖性）
- 5-HTP 是 5-HT 前体，临床上用于抑郁/纤维肌痛（Maffei 2020, IJMS 综述）
- BBB 5-HTP 转运机制**不明确**：LAT1 介导（Kanai 1998）+ 内皮脱羧酶捕获（Hardebo 1979）被提出，但无专属转运体确证
- 抑郁症患者 [11C]5-HTP 脑摄取降低（Agren 1991, PET 研究）提示转运缺陷
- **SMVT 作为 5-HTP 候选转运体 = 本稿件假说的新颖点**，尚无直接实验证据 → 需实验验证（5-HTP 竞争摄取实验）

## 稿件引用状态
- `uchida2015` — **DOI 已修正** (10.1111/jnc.13089 → 10.1111/jnc.13092)
- `shaikh2023` — 验证准确（累积 NSAID → 抑郁正相关）
- 待补：`hu2020` (BMC Med, 非阿司匹林 NSAID HR 1.24) 可选作更强流行病学证据

## 待验证实验（稿件 Limitations 已列）
1. 5-HTP 竞争摄取实验（SMVT 表达细胞系）
2. NSAID 对 SMVT 介导 5-HTP 摄取的抑制实验
3. BBB 模型（hCMEC/D3 或小鼠）验证 NSAID→5-HTP→5-HT 通量
