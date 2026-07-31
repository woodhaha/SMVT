# Research Log — SMVT

## 2026-07-31 — Autoresearch 启动（稿件 figure 精修阶段）
- 状态：已有完整 8 化合物 MD + 泛癌 DE + 突变景观。进入稿件 figure 精修。
- 进行中工作：Fig1 (pan-cancer DE) 重写完成跑通；Fig3 (mutation landscape) 在改。
- **问题**：fig_mutation_final.py 读 `mc3_maf.gz`（不存在），只有未压缩 3.5G `mc3.v0.2.8.PUBLIC.maf`。
- **修复**：写 `extract_slc5a6_maf.py` 预提取 SLC5A6 行 → `02_Data/cleaned/slc5a6_mc3_maf.csv`（后台跑）；改脚本读 CSV。
- 决策：预提取后 mutation 图迭代从 ~2min/次 降到 ~1s/次。

## 待做
- [ ] 验证 slc5a6_mc3_maf.csv 提取完成 → 重跑 fig_mutation_final.py
- [ ] 检查 Fig1/Fig3 输出 → outer loop 排下一个缺口
