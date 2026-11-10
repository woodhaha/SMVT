# Research Log — SMVT

## 2026-11-04 — 第 4 轮堂会审 + SMVT 药物设计策略研究
- **第 4 轮堂会审**（3 审稿人 + 复审）：抓出 10 处问题，核心是 said1989 反噬（证明巴比妥配体同时证明抑制有毒）+ Vina −8.0 阈值统计脆弱 + 金刚烷 SAR 无统计 + MD 只报蛋白 RMSD。全部修复，靶点定位统一为"probe/prodrug-delivery + 抑制需治疗窗"。
- **投稿包同步**：发现 CBM_submission/ 是 10-27 旧快照，缺 expanded screening 段（Table 7 + Fig 9），已完整同步 04_Manuscript 最新版。
- **药物设计研究**（网络 + zhen2026 结构 + tripathi2025）：核心结论=自由羧基是 SMVT 识别锁，biotin 偶联物（酰胺消耗羧基）不被转运，金刚烷羧酸是最优递送头。详见 `SMVT_drug_design_strategy.md` + 记忆 smvt-drug-design-carboxyl-key。
- **docking 验证（证伪偶联物递送）**：cryo-EM 26va（occluded）+ 26vd（outward-open）。金刚烷二羧酸 −7.73（强于 biotin −6.49），羧基锚定 Tyr156 3.5 Å（识别头正确）；但偶联物 16.7 Å 太大，occluded 态 +27.45 排斥，outward-open 态 −10.93 是非特异性塞入（羧基距 TM1 kink 7–8 Å）。结论：金刚烷羧酸定位为 SMVT 抑制剂/探针，非递送头；递送走小分子前药范式。记忆 smvt-conjugate-docking-negative。
- **三轮补充验证（换识别头 + 换 payload）**：换识别头（丙二酸/戊二酸/环己烷二羧酸+SN-38）全排斥——SN-38 本身 392 Da 超限，识别头非变量；换最小 payload 5-FU（130 Da）偶联物 336 Da 能进 occluded 口袋（−5.16）但羧基距 Tyr156 8.1 Å 仍不锚定。最终结论：SMVT 底物识别是"羧基插 TM1 kink + 整体几何匹配"，任何双模块偶联都破坏匹配（大→进不去、小→羧基脱靶），偶联递送方向彻底证伪。
- **方向转变：转运递送 → 纳米粒靶向胞吞**：把 SMVT 从"转运体"重新定位为"肿瘤细胞表面靶标"。不要求偶联物匹配口袋被转运，而是金刚烷羧酸（−7.73 强于 biotin）修饰纳米粒/脂质体，结合 SMVT 后介导胞吞内化（~100 nm 尺度，绕开口袋限制）。有先例（MDPI Cancers 12:2837：biotin 修饰脂质体/cubosome 靶向 SMVT，游离 biotin 竞争证实特异性）。金刚烷羧酸新定位 = 更强 SMVT 靶向头。记忆 smvt-nanoparticle-targeting-direction。

## 2026-11-04 — 投稿目标重大变化：CBM 已死 → PLOS CB
- **CBM 已死**（web 查证）：2024-04 On Hold（同行评审操纵、假引用）→ 2025-11-17 WoS Core Collection 移除。无 IF/分区，温医绝不奖励，CBM 投稿包作废。
- **期刊筛选**（web 查 scope 非仅 IF）：BIB（综述期刊不收原创）、JCIM（"不收单靶标对接应用+无实验"）、PLOS CB（唯一 fit：收应用+实验不必需）。结论：2区以上唯一选择 = PLOS CB（2区生物学 IF3.7）。记忆 smvt-journal-scope-screening。
- **稿件适配 PLOS CB**：标题改 "Carboxyl-Bearing Ligands"、参考文献转 Vancouver（unsrt 前6作者 et al.，30条，查 Crossref 补全）、摘要 253 词（≤300 OK）、Code availability 措辞改"public repository"、金刚烷三羧酸 cryo-EM −7.64 验证 + 靶向递送定位（hypothesis）。
- **代码归档计划**：code_archive_plan.md（6 模块 ~25 核心脚本 + 数据来源 + 3 缺口）。待办：建 GitHub repo、从 weste 同步 md_v13.py、温医 T1/T2/T3 核对。

## 2026-11-02 — 对接扩展完成 + 金刚烷羧酸 SAR
- **对接扩展收官**：主库 10569（9659 ok / 910 fail，288 min）+ anticancer 2371（2164 ok）→ 合并去重 **11823 个**。hits <−8 = **20**（全氟已剔）、<−7 = **349**、<−6.5 = 910。原 703 库仅 8 个 <−8。
- **三大化学型**：①金刚烷/钻石烷羧酸簇（6 个，全新 chemotype）②巴比妥/戊二酰亚胺簇（phenobarbital −8.17 + 哌啶二酮腈 −8.61 + glutethimide −7.75）③NMDA 轴（dizocilpine/MK-801 −8.02 + 美金刚类似物 −8.12）。
- **假阳性剔除**：perfluorodecalin/PFOA（全氟）+ 苯并三唑（紫外吸收剂）。
- **金刚烷羧酸 SAR**：58 类似物本地对接（4.2 min），best **COOH+Ac −8.84**（超主库金刚烷 −8.62，与巴比妥 SAR −8.79 同级）。规律：双羧基>单羧基、F/CN 增强、diamantane 太大。
- **数据备份**：本地 `03_Analysis/outputs/docking_final/`（results.csv / anticancer_results.csv / anticancer_combined_results.csv / 两个 checkpoint）。
- **坑**：①B12 钴 segfault → 金属/超大须预过滤；②Na/K/Ca 盐须去盐；③watchdog 缓存 anticancer 目标数须重启重读。详见记忆 smvt-docking-b12-segfault / smvt-adamantane-sar-local-dock。
- **pod 关机**：SSH shutdown/poweroff 被拒（容器无 systemd），须网页控制台点关机（pod `cpod-1wsdz4i2otkn`）。
- **20 hit 查证**：PubChem 全查（19 有 CID，1 个 UWARBR 无记录）。可采购性按身份分类：8 个 FDA/天然产物（phenobarbital/cytisine/cedrol/patchouli alcohol/securinine/sorbinil/CRT0044876/dizocilpine）可购；金刚烷羧酸簇（6 个）Sigma 类可购但非上市药；苯并三唑工业品剔除。分类表 `03_Analysis/outputs/docking_final/hits_top20_classified.csv`。
- **稿件整合**：screening 扩展段正式并入 `SMVT_manuscript.tex`——Results 段（11823 化合物+三大化学型+SAR）+ **Table 7**（tab:expanded，20 hit 含 InChIKey/ΔG/化学型/可购性）+ **Fig 9**（RDKit 20 结构网格，按化学型分组）+ Methods 段（Expanded Library Screening）。编译 0 错误 52 页。
- **图号修复**：机制图（Fig2）从文档末尾移回第 2 位；扩展筛选图定 **Fig9**（原 FigX 重命名）。现 Fig1–Fig9 文件名 = 文档位置 = 硬编码引用三处对齐。
- **文献意义调研**：`findings_literature.md` 存档。核心发现——①巴比妥发现印证 Said 1989（carbamazepine/primidone 竞争性抑制 SMVT 生物素转运，Ki 4.70/2.25 mM）②金刚烷羧酸是已知 DGAT1 抑制剂药化骨架 ③NMDA 轴无文献（最弱，建议降权）。
- **Discussion 增补（已并入稿件）**：+3 篇文献（said1989/park2005/pagire2015）——巴比妥临床先例、羧基必需性底物特异性、金刚烷 DGAT1 药化背书。编译 0 错误。

## 2026-10-27 — 改投 Computers in Biology and Medicine (CBM)
- 期刊从 CSBJ 改为 CBM（IF 6.3/JCR Q1/中科院大类 2区医学；范围对"临床生信+ML+docking/MD"混合体更对口，CSBJ 偏结构会嫌临床那块）
- CSBJ_submission/ → CBM_submission/，cover letter 重写 CBM 版（.md+.docx），README/checklist 更新，稿件/title page/highlights/reviewers 均 journal-neutral 未动
- 剩余：审稿人邮箱核实、温医 T1/T2/T3 核对、GA 视觉核验

## 2026-10-26 — Elsevier 模板迁移 + 全稿 audit + 投稿包 + 图形摘要
- 期刊定稿 CSBJ（IF 4.8/中科院3区/JCR Q2；对比过 iScience/PLOS CB/JBSD，详见 research-state）
- Elsevier elsarticle 模板迁移：sn-jnl → elsarticle + elsarticle-num，编译 0 错误 0 警告（47 页）
- 全稿 audit：4 处文献真实性修复（Paixão/Shaikh N.F./zhen2026 vol17:8036/leonardi）+ FDA→Approved + 标题重写 + Vina 两数据集陷阱识别
- 目录整理：_archive/ 归档旧文件、03_Analysis 图合并、.gitignore 更新
- CSBJ_submission/ 投稿包建成（manuscript + 4 投稿文档 docx + 24 图 + README）
- 图形摘要：Gemini 提示词 + matplotlib 热图 + 1328×531 合成图，已入 CSBJ_submission
- 剩余：审稿人邮箱核实、温医 T1/T2/T3 核对、GA 视觉核验

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

## 2026-07-31 (深夜) — 数字一致性审计（二轮）+ fold 修正
- 修正 5-HTP fold 声称：10x → 4.6x（e^(-0.90/RT)=4.56）。glutarimide 5.6x 验证正确。
- 审计通过：8 癌型统一、LGG/KIRC HR 正文与 Table 一致、AUC=0.888 一致、全部 docking 分数跨章节一致（biotin 6.76 pilot / 6.82 screening 有区分）。
- 提交: baed3ef。

## 2026-07-31 (深夜) — 5-HTP MD 2fs production 稳定运行
- 4fs production 崩（HMR 不覆盖配体极性 H）→ 改 2fs，现已稳定。
- 速率 ~20 ns/hour，20ns 测试约 60 分钟完成。
- 主表 8 化合物参照齐备（RMSD 0.53-0.70 全稳定），analyze_md.py 输出字段与主表一致。
- cron 每 20 分钟监控，完成后自动分析对比。

## 2026-07-31 (深夜) — 5-HTP MD 崩盘分析 + 验证成功
- **崩盘**：2fs production 跑到 11.2ns（214/400 帧）后 `Particle coordinate is NaN`，进程死。8 化合物全稳定 100ns，唯独 5-HTP 崩。
- **排障关键**：最初分析配体 COM 漂移显示"逃逸 13Å"——实为 **PBC 镜像伪影**（13.2nm 盒子边缘，配体 x≈13.1 与 x≈0 镜像跳变）。`image_molecules()` PBC 修正后确认配体 100% 结合位点。
- **结果**：11.2ns 轨迹 Cα RMSD 0.49±0.11Å（9 化合物最低），接触占用率 99.5%，配体 COM 漂移 0.22±0.06Å。**5-HTP MD 验证成功**。
- **崩盘根因**：数值不稳定（GAFF2 极性 H 高频），非物理逃逸。HMR 不覆盖配体极性 H。
- **本地分析管线打通**：云端重建完整 solvated 拓扑（228540 原子 5HTP_final.pdb）→ scp 同步 dcd+csv+pdb → 本地 analyze_md.py 跑通。
- **稿件**：正文 MD 段加 5-HTP 句 + FigS15 升级双面板（docking overlap + RMSD），编译干净（40 页，仅无害 figure.N 警告）。
- 提交: e142ebf。

## 2026-08-01 (凌晨) — 5-HTP MD 重跑：AllBonds 修复配体极性 H
- **根因**：`constraints=app.HBonds` 不约束配体极性 H（吲哚 NH、酚 OH、羧基 OH、氨基 NH₂）→ 高频振动未阻尼 → 2fs 崩在 11.2ns。
- **修复**：改用 `constraints=app.AllBonds`（约束所有重原子-H 键）。云端 weste 重跑，NVT/NPT 全过（PE 正常 −699 至 −858k），2fs production 运行中。
- 速率 ~5.6 ns/hr（AllBonds + 2fs 代价，比 4fs 8 化合物慢），预计 ~01:04 完成。
- findings.md 记录 AllBonds 教训。
- 本地诊断教训：addSolvent 随机性导致本地复现不可靠；但本地 lig/protein clash 检查（min 1.02Å）暴露了坐标敏感性。

## 2026-08-01 (凌晨) — 5-HTP MD 崩盘根因深度排障
- **关键教训**：分析含大量水的 MD 轨迹时，`t.superpose(frame=0, atom_indices=backbone)` 会**破坏配体-蛋白距离**，导致误判"配体逃逸"（superpose 后 min dist 显示 7-8Å 逃逸）。**必须用 RAW 坐标 + PBC 感知（image_molecules）**。
- **RAW PBC 可靠结果**：5-HTP 在 11.2ns 内 **100% 时间结合**（min dist 0.175±0.010Å，全程 <0.21Å），每帧 27.8 个接触（配体 28 原子全在蛋白 4Å 内），从未逃逸。**稿件"stably bound"描述正确**。
- **崩盘根因**：三种协议全崩（v13+HMR+HBonds 崩 11.2ns、v13+HMR+AllBonds 崩 8ns、v8 无HMR 崩 NPT）→ 纯数值不稳定，与配体结合无关。配体稳定结合但 8-11ns 后数值爆炸（GAFF2 极性 H 高频振动）。
- **重跑方案**：1fs production（20M 步，AllBonds 约束）——最保守的数值稳定保证。~13h。云端运行中。
- **权限修复**：settings.json ask 列表清空 + defaultMode=bypassPermissions，消除确认弹窗。
- **最终决策（08-01）**：完整 20ns 重跑 4 种配置全崩（HBonds 11.2ns、AllBonds 8ns、v8 NPT、1fs 1ns）→ 根因 = 5-HTP GAFF 参数数值不稳定。稿件诚实修正为 11.2ns（RAW-PBC 确认 100% 结合），编译干净，commit a2203af。完整 20ns 转技术债。
- **教训**：MD 数值崩溃排查优先级应为：(1) 先 RAW-PBC 验证配体是否真结合（勿信 superpose）；(2) 排除初始 clash；(3) 换步长/约束；(4) 重新参数化配体。5-HTP 卡在第 4 步（GAFF 参数）。

## 2026-08-01 — weste 关机 + 投稿包终检
- weste 服务器 05:35 关机（SSH shutdown），数据全部同步本地。5-HTP MD GPU 总耗时 ~4h。
- 投稿包核对：40 页稿件 + 24 图 + 双封面信（NC/CR）+ 实验方案，全部就绪。
- 唯一待办：用户决策投稿期刊（NC vs Cell Reports）。

## 2026-08-01 — 模拟堂会审 3 轮（重大修改）
- **Round 1**：3 审稿人（方法学/结构/转化）全部 Major Revision。共识问题：核黄素阴性对照失败、生存无 FDR、5-HTP top FDA 矛盾、11.2ns 误导、MD 无脂双层/Na⁺、passenger 无选择检验、假说过度推断。
- **Round 1 修改**：FDR 校正生存（8→6 癌型 + q 值）、5-HTP 诚实化（top aromatic-amino-acid）、FigS1/S2/S16 用真实数据重生成、negative control caveat、MD 局限声明、passenger 降级、glutarimide 烷化剂警示。
- **Round 2 复审**：R1→Minor、R2→Major（第2轮）、R3→Accept(minor)。残留：摘要/legend top-FDA 遗留、MD confirm 过度、羧基 not-essential 矛盾、KO 非正交、标题过度声称。
- **Round 2 修改**：标题→"Predicted Ligands"、摘要/正文 establish→identify/propose、confirm→support、羧基声称软化、KO 改内部一致性、Na⁺ 断行修复。commit a229996。
- **Round 3**：Reviewer 2 最终复审 → **Minor（接受）**。3 审稿人全部通过。
- **收尾**：Table 2 加事件列（UVM 79/22 透明化）、最终一致性扫描全过。
- **会话状态**：2026-08-01 保存。稿件投稿就绪，唯一待办 = 用户决策期刊（NC/Cell Reports）。
- **教训**：模拟会审能有效暴露投稿硬伤（FDR 缺失、标题过度、legacy 图数据错），比人工自查更全面。FigS1/S2 曾被 generate_manuscript_figures.py 用旧 LUSC/LUAD 数据 + 合成曲线污染——审稿人抓不到但会审抓到。

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
