# SMVT MD — Quick Reference

> Last updated: 2026-07-01 00:12 CST · 8 compounds · RTX 5090

## Connect
```bash
ssh -o ProxyCommand=none -p 39591 root@connect.bjb2.seetacloud.com
```

## Monitor
```bash
# GPU
ssh -p 39591 root@connect.bjb2.seetacloud.com nvidia-smi

# All logs
ssh -p 39591 root@connect.bjb2.seetacloud.com "tail -f /root/autodl-tmp/SMVT_AutoDL_Package/logs/*.log"

# Single compound
ssh -p 39591 root@connect.bjb2.seetacloud.com "tail -20 /root/autodl-tmp/SMVT_AutoDL_Package/logs/BIOTIN_v13.log"
```

## Files
| Path | Content |
|------|---------|
| `/root/autodl-tmp/SMVT_AutoDL_Package/` | Data root |
| `scripts/md_v13.py` | MD engine (HMR+4fs) |
| `run_parallel_v2.sh` | Parallel launcher |
| `logs/{NAME}_v13.log` | Per-compound log |
| `trajectories/{NAME}/` | DCD+CSV output |
| `ligands/{NAME}.pdbqt` | Vina docking pose |
| `ligands/{NAME}_template.xml` | GAFF forcefield template |
| `receptor/AF-Q9Y289-F1.pdb` | SMVT protein |

## Environment
- **Conda**: `md_env` (Python 3.11.15)
- **OpenMM**: 8.5.2 source-built, CUDA 12.8, sm_120
- **Activate**: `source /root/miniconda3/etc/profile.d/conda.sh && conda activate md_env`

## 8 Compounds
| # | Name | Role |
|---|------|------|
| 1 | NAFTAZONE | Test (best Phase A) |
| 2 | RIBOFLAVIN | Negative control (vitamin B2) |
| 3 | FUROSEMIDE | Test (WHO essential) |
| 4 | PHENOBARBITAL | Test (barbiturate) |
| 5 | HYDROMORPHONE | Test (best hit, dG -8.58) |
| 6 | GABAPENTIN_ENACARBIL | Positive control (SMVT prodrug) |
| 7 | ESKETAMINE | Test (NMDA antagonist) |
| 8 | BIOTIN | Reference (natural substrate) |

## Protocol (v13)
1. PDBFixer → 2. Vina pose → 3. RDKit ligand → 4. ForceField(amber14+GAFF2+TIP3P) → 5. Solvate → 6. HMR → 7. Backbone restraints → 8. Minimize 5K → 9. NVT 50→300K (2fs) → 10. NPT equil (2fs) → 11. Production 100ns (4fs)

## Local Files
- `D:\Researching\SMVT\SMVT_AutoDL_Package\` — full package copy
- `C:\Users\woodh\Documents\md_v13.py` — MD script
- `C:\Users\woodh\Documents\run_parallel_v2.sh` — launcher
