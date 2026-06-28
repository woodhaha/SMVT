# [蒸馏] scTenifoldKnk 虚拟敲除 + 孟德尔随机化 联合分析范式

> Source: 科研解忧铺 · 云生信-百合 · 2026-06-20
> URL: https://mp.weixin.qq.com/s/35MTFUeiaROjHK463zcjuQ
> Distilled: 2026-06-23 · 对 SMVT 项目有直接方法论参考价值

---

## 1. 文章概要

| 维度 | 详情 |
|------|------|
| **标题** | 9成生信！单细胞+机器学习哪够，再来scTenifoldKnk虚拟敲除+孟德尔随机化做验证 |
| **团队** | 安徽中医药大学 |
| **靶基因** | **CKAP2** (cytoskeleton-associated protein 2) |
| **疾病** | 胃癌 (GC) |
| **投稿→接收** | **3 个月** |
| **生信占比** | ~90%（仅最后做 WB 实验验证） |

---

## 2. 完整分析管线（10 步）

```
Step 1: DEGs + WGCNA
  ├── 训练集: TCGA-STAD
  ├── 验证集: GSE13911, GSE65801, GSE13861, GSE29272
  └── 取交集 → 候选基因集

Step 2: 12 种机器学习算法 → 113 个集成模型
  └── 按 C-index 筛选最佳模型 → 筛选候选 genes

Step 3: STRING PPI + 网络拓扑分析
  └── 核心枢纽 genes

Step 4: 单细胞 RNA-seq (GSE163558)
  ├── 标准化 → 聚类 → 注释细胞类型
  └── 确定关键基因在细胞内的定位

Step 5: 上皮细胞亚聚类 + inferCNV
  └── 鉴定恶性上皮细胞

Step 6: ⭐ scTenifoldKnk 虚拟敲除
  ├── 构建基因调控网络 (GRN)
  ├── 扰动 CKAP2 → 比较扰动前后网络
  ├── 识别差异调控基因 (DRGs)
  └── DRGs 富集分析 → 5-HT 通路

Step 7: CellChat 细胞间通讯分析

Step 8: 临床关联分析
  ├── Wilcoxon/KW: CKAP2 vs 临床病理特征
  └── Spearman: CKAP2 vs TMB

Step 9: ⭐ 双向孟德尔随机化 (MR)
  └── 评估 CKAP2 与 GC 的因果关系

Step 10: WB 实验验证
  └── CKAP2 在 GC 组织和细胞系中上调
```

---

## 3. 关键发现

| 发现 | 方法 |
|------|------|
| CKAP2 在 GC 恶性上皮细胞中高表达 | scRNA-seq + inferCNV |
| CKAP2 是独立保护性预后标志物 | ML + 生存分析 |
| **虚拟 KO CKAP2 → DRGs 富集于 5-HT 通路** | scTenifoldKnk |
| CKAP2 与 GC 存在因果关系（保护方向） | 双向 MR |
| CKAP2 在 GC 组织和细胞系中上调 | WB |

---

## 4. scTenifoldKnk 方法要点

**scTenifoldKnk** = 单细胞基因调控网络 (GRN) 虚拟敲除工具

| 步骤 | 操作 |
|------|------|
| 1 | 从 scRNA 数据构建基因调控网络 (GRN) |
| 2 | 从网络中"移除"目标基因（虚拟 KO） |
| 3 | 比较 KO 前后网络拓扑变化 |
| 4 | 识别差异调控基因 (Differentially Regulated Genes, DRGs) |
| 5 | 对 DRGs 进行通路富集 → 推断 KO 的功能后果 |

**优势**：无需湿实验 CRIPSR，利用已有 scRNA 数据即可预测基因 KO 的转录组后果。

---

## 5. 对 SMVT 项目的可迁移性

| 文章 CKAP2 管线 | SMVT 可迁移性 | 状态 |
|-----------------|-------------|:---:|
| DEGs + WGCNA | ✅ 已有 TCGA 表达数据 | 已完成 |
| 12 ML 算法筛选 | ⚠️ 可补充（当前仅有表达分析，无 ML 模型） | 待做 |
| scRNA-seq 定位 | ⚠️ 需获取 GC 或其他癌种 scRNA 数据 | 待做 |
| **scTenifoldKnk 虚拟 KO** | ✅ **已完成虚拟 KO 推理分析**，可升级为定量 GRN-KO | 可升级 |
| CellChat 细胞通讯 | ⚠️ 需 scRNA 数据 | 待做 |
| 孟德尔随机化 | ⚠️ 可补充 SMVT→癌症因果关系 | 待做 |
| WB 实验验证 | ⚠️ 需湿实验资源 | 外部 |

---

## 6. 数据资源速查

| 资源 | ID |
|------|-----|
| TCGA 胃癌训练集 | TCGA-STAD |
| 验证集 1 | GSE13911 |
| 验证集 2 | GSE65801 |
| 验证集 3 | GSE13861 |
| 验证集 4 | GSE29272 |
| 单细胞数据 | GSE163558 |
| 虚拟敲除工具 | scTenifoldKnk (R package) |

---

## 7. 核心 Takeaway

> **"DEGs + WGCNA + 12 ML + scRNA + scTenifoldKnk 虚拟 KO + MR + WB" = 90% 生信 + 10% 实验 = 3 个月投稿→接收**
>
> 这个范式对 SMVT 项目高度可迁移。当前 SMVT 已完成 DEGs/WGCNA 对应步骤（表达 + 突变 + STRING + 富集），缺失的是：(1) ML 模型筛选、(2) scRNA 单细胞定位、(3) 定量 GRN 虚拟 KO、(4) MR 因果验证。
>
> **scTenifoldKnk 是当前最直接的可升级项** — 将 SMVT 虚拟 KO 从定性推理升级为定量网络扰动分析。
