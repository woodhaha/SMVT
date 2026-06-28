# SLC5A6 (SMVT) 结构准备报告

> 日期: 2026-06-23 · 工具: AlphaFold DB, Biopython, OpenBabel, OpenMM

---

## 1. 结构来源

| 项目 | 详情 |
|------|------|
| **蛋白** | Sodium-dependent multivitamin transporter (SMVT) |
| **基因** | *SLC5A6* |
| **UniProt** | Q9Y289 (SC5A6_HUMAN) |
| **来源** | AlphaFold v6 (无 PDB 实验结构——膜蛋白结晶困难) |
| **长度** | 635 氨基酸 |
| **全局 pLDDT** | 79.38 |
| **置信度分布** | VH(>90): 51.2%, Confident(70-90): 27.6%, Low(50-70): 3.9%, VL(<50): 17.3% |
| **下载 URL** | `https://alphafold.ebi.ac.uk/files/AF-Q9Y289-F1-model_v6.pdb` |

---

## 2. 结构准备流程

### Step 1 — 下载
```
AF-Q9Y289-F1.pdb  (398,681 bytes, 4,824 atoms, 635 residues)
```

### Step 2 — 去除水分子和杂原子
```
工具: Biopython PDBParser + Select
结果: 0 个水分子, 0 个杂原子 (AlphaFold 预测结构本身不含)
输出: AF-Q9Y289-F1_cleaned.pdb (395,659 bytes)
```

### Step 3 — 加氢原子 (pH 7.0)
```
工具: OpenMM Modeller
方法: modeller.addHydrogens(forcefield, pH=7.0)
力场: Amber14 protein.ff14SB
结果: 4,824 → 9,746 atoms (+4,922 H)
输出: AF-Q9Y289-F1_H.pdb (~1.5 MB)
```

### Step 4 — 力场电荷
```
工具: OpenMM (Amber14SB 内置)
方法: Amber ff14SB 部分电荷 (力场参数内嵌, 非 Gasteiger)
注: Amber14SB 电荷比 Gasteiger 更适合蛋白质模拟
```

### Step 5 — 能量最小化 ✅
```
工具: OpenMM 8.5.2
力场: Amber14SB + GBSA-OBC 隐式溶剂
方法:
  Phase A: 重原子约束 100 kcal/mol/Å² → H 弛豫 (300 步)
  Phase B: 全部释放 → 全原子最小化 (2000 步)
结果:
  初始能量:  -39,870  kJ/mol
  H-弛豫后:  -51,667  kJ/mol
  最终能量:  -55,275  kJ/mol
  总降幅:    15,405  kJ/mol (3.19 kJ/mol per 重原子)
输出: AF-Q9Y289-F1_prepared.pdb (~1.5 MB)
```

---

## 3. 输出文件清单

| 文件 | 大小 | 状态 |
|------|------|------|
| `AF-Q9Y289-F1.pdb` | 399 KB | 原始 AlphaFold v6 |
| `AF-Q9Y289-F1_cleaned.pdb` | 396 KB | 去杂后 (无 H) |
| `AF-Q9Y289-F1_H.pdb` | 1.5 MB | OpenMM 加氢 (pH 7.0) ✅ |
| `AF-Q9Y289-F1_prepared.pdb` | 1.5 MB | 能量最小化 ✅ |

---

## 4. SMVT 关键残基

基于文献的功能残基（用于对接盒子定义）：

| 功能 | 残基 |
|------|------|
| Na+ 结合位点 1 | Ser-?, Asn-?, Thr-? (待从同源建模确认) |
| Na+ 结合位点 2 | 同上 |
| Biotin 结合口袋 | 跨膜螺旋 TM1-TM6 的中央腔 |
| 保守 motif | SLC5 家族 Na+/溶质 symport 共有序列 |

> ⚠️ 注: AlphaFold 预测结构中确切的底物结合位点需与实验结构 (如 vSGLT 同源模板，PDB: 3DH4) 做结构比对确认

---

## 5. 下一步

1. [x] ~~OpenMM 能量最小化 (Amber14SB + GBSA 隐式溶剂)~~ ✅
2. [ ] 与 vSGLT (3DH4) 结构对齐 → 定位底物结合口袋
3. [ ] 对接盒子定义 + 分子对接 (AutoDock Vina / Smina)
4. [ ] MD 模拟 (GROMACS / OpenMM) 验证稳定性
