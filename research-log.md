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
