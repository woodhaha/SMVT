# Research Log — SMVT

## 2026-07-31 — Autoresearch 启动（稿件 figure 精修阶段）
- 状态：已有完整 8 化合物 MD + 泛癌 DE + 突变景观。进入稿件 figure 精修。
- 进行中工作：Fig1 (pan-cancer DE) 重写完成跑通；Fig3 (mutation landscape) 在改。
- **问题**：fig_mutation_final.py 读 `mc3_maf.gz`（不存在），只有未压缩 3.5G `mc3.v0.2.8.PUBLIC.maf`。
- **修复**：写 `extract_slc5a6_maf.py` 预提取 SLC5A6 行 → `02_Data/cleaned/slc5a6_mc3_maf.csv`（后台跑）；改脚本读 CSV。
- 决策：预提取后 mutation 图迭代从 ~2min/次 降到 ~1s/次。

## 待做
- [x] 验证 slc5a6_mc3_maf.csv 提取完成 → 重跑 fig_mutation_final.py
- [x] 检查 Fig1/Fig3 输出 → outer loop 排下一个缺口

## 2026-07-31 — 稿件一致性审计完成
- 逐数字核对稿件 vs 数据，修两处不一致：DE (Table1: 6→14癌型, log2FC改新口径, 统计方法改Mann-Whitney U) + mutation (141/82→193/114, V120I 4→5, burden 0.6%→1.5%)。
- 其余数字全部验证正确：对接分、MD、ML AUC。
- 修 LaTeX 预存错误（缺 amsmath），稿件编译干净。
- 提交: c2cf1d8 fix(manuscript)。

## 下一步候选
- 突变致病性 (ESM-1v/αMissense) — 稿件已提方法但输出未整合?
- MD↔对接 ΔG 一致性 (Furosemide 偏差)
- 泛癌生存合成进叙事

## 2026-07-31 (晚) — FigS1/S2 生存图一致性修复
- 发现 FigS1（4 癌型 LUSC/LUAD/BLCA/LIHC）和 FigS2（10 癌型）用的是旧 cBioPortal 数据，与 Table 2（8 癌型 LGG/KIRC/ACC/MESO/KIRP/SARC/UVM/STAD）矛盾。
- 从 `pancan_survival.csv`（已验证匹配 Table 2）重生成：FigS1 = 8 面板复合 KM + HR/P 标注；FigS2 = 8 癌型 forest plot（log scale）。
- caption 更新 4/10 → 8 癌型。编译干净。
- 提交: 6b2d02a, 564257e。
- FigS3 验证：数据源（84 化合物 docking_expanded_results.csv）未变，图正确。

## 2026-07-31 (深夜) — HPA 蛋白验证：SLC5A6 肾癌不利预后标志物
- HPA API + WebSearch 确认：SLC5A6 蛋白水平证据 + 质膜定位（surface drug target 支撑）。
- **HPA 明确列 SLC5A6 为肾癌不利预后标志物**——独立于 TCGA 的蛋白层验证，与稿件 KIRC HR=9.2 吻合。
- 已加入稿件 122 行生存段 + uhlens2015 引用（DOI 验证）。编译后台进行。

## 2026-07-31 (深夜) — SMVT 结合特异性：5-HTP vs SSRI 对照
- 云端 GPU 未就绪，5-HTP MD 等待中 → 转做零 GPU 的稿件补强。
- 发现稿件血清素轴假说缺 SSRI 对照组，审稿人会问"SMVT 是否非特异"。
- 同一 R1-R3 轮数据：5-HTP −7.66 >> Fluoxetine −6.13 > Sertraline −4.25。SMVT 偏好 5-HTP 而非 SERT 抑制剂本身。
- 已加入稿件 Discussion。编译后台验证中。

## 2026-07-31 (夜) — ESM-2 650M 交叉验证（passenger 结论稳健性）
- 动机：passenger 结论只靠单一 35M 模型，审稿人可挑战。
- 用 HF 缓存的 esm2_t33_650M_UR50D 重算 87 变异（GPU，秒级）。
- **结果**：35M vs 650M ρ=0.806 (p=5e-21)；R123L 对照 650M 排名 13%（更强致病）；650M 最致病 18 变异 17 个 >10Å 离 biotin 位点；分布无偏无簇集。
- 稿件 Results/Methods 加 robustness 描述，编译干净。
- 提交: (本轮)

## 2026-07-31 (深夜) — MC3 barcode map 完整化 + Fig3/Fig4 重生成
- **关键 bug 修复**：code_map 漏 18 个 TCGA 中心代码，丢 22 个突变。补全后 Fig3 显示 193/114/167（完整）。
- **正文修正**：mutation 频率 "<2% all" → 准确表述（ESCA 5.9% 等超 2%，随高突变癌型基线）。
- Fig3/Fig4 用唯一患者/标准 TCGA 样本数重生成，caption 更新。
- 提交: c0757fb。
