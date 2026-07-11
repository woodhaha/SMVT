# SMVT MD — Research Plan & Timeline

> Created: 2026-06-29 | Updated: 2026-06-29

## Current Status

| Phase | Status | Detail |
|-------|--------|--------|
| Phase 1: MD Production | 🟢 Running | 8 compounds × 100ns, RTX 5090 |
| BIOTIN | 🟢 7.3/100ns done | ~10.4h remaining (v9, 2fs, 214 ns/day) |
| PHENOBARBITAL → ESKETAMINE (7) | 🟡 Queued | v11 (4fs HMR + mixed, ~515 ns/day) |
| Estimated completion | — | ~43h total (~1.8 days) |

## Technical Stack

| Component | Detail |
|-----------|--------|
| GPU | NVIDIA RTX 5090, 32GB, CC 12.0, CUDA 13.0 |
| MD Engine | OpenMM 8.5.2 (source-built, CUDA 12.8, sm_120) |
| Force Field | AMBER14 + GAFF2 + TIP3P |
| Protein | SMVT (AF-Q9Y289-F1) |
| Environment | conda md_env, Python 3.11 |
| Server | SeetaCloud (connect.westd.seetacloud.com:16492) |

## Optimization History

| Version | Changes | ns/day | xBase |
|---------|---------|--------|-------|
| v8 | Base (2fs) | ~52 | 1× |
| v9 | PDBQT fix + barostat always-on + NPT fix | ~214 | 4× |
| v10 | HMR (H mass 3amu) + 4fs timestep | ~430 | 8× |
| v11 | Mixed CUDA precision | ~515 | 10× |

Key fixes:
- PDBQT parser: supports ROOT/BRANCH format (was MODEL/ENDMDL only)
- NPT: barostat in system from start, no context rebuild (reinitialize + addForce caused NaN)
- HMR: adjacency-list O(n+m) lookups, hydrogen mass 1→3 amu

## Publication Strategy

### Phase-by-Phase

#### Phase A: Pure MD (Current — Low-mid tier)
- 8 compounds × 50ns MD + RMSD/RMSF/Rg analysis
- Target: J Mol Model (IF 2-3) / Comput Biol Chem (IF 1-3)
- Minimum publishable unit

#### Phase B: MD + Binding Free Energy (Mid tier)
- Add MM-GBSA/PBSA per-residue decomposition
- Mutation hotspot prediction (alanine scanning in silico)
- Target: JCIM (IF 4-5) / J Chem Inf Model
- **Recommended next step after Phase 1**

#### Phase C: MD + Experiment (Mid-upper tier)
- SPR/ITC/MST binding affinity validation
- Mutation + functional assay
- Target: J Med Chem (IF 6-8) / Eur J Med Chem

#### Phase D: Full Story (High tier)
- Cryo-EM or NMR structure
- Functional experiments
- Target: Nature Commun / Sci Adv / eLife (IF 10+)

### Story Angles

1. **Old protein, new ligands**: SMVT is a known transporter; drug repurposing
2. **BIOTIN as natural substrate control** vs 7 drug candidates
3. **Computational-to-experimental pipeline**: MD → MM-GBSA → experimental validation

## Analysis Pipeline (After MD completes)

```bash
# 1. RMSD/RMSF
python analysis/rmsd_rmsf.py trajectories/*/

# 2. MM-GBSA binding free energy  
python analysis/mmgbsa.py trajectories/*/

# 3. Contact maps / hydrogen bond analysis
python analysis/contacts.py trajectories/*/

# 4. PCA / clustering of conformations
python analysis/pca_cluster.py trajectories/*/
```

## Next Actions

1. ✅ BIOTIN complete → auto-trigger 7 compounds
2. ⬜ Verify all 8 compounds finish without NaN
3. ⬜ Download trajectories + CSVs locally
4. ⬜ Run MM-GBSA analysis
5. ⬜ Draft manuscript abstract + figures
6. ⬜ Submit to JCIM

## Key Commands

```bash
# Check progress
ssh autodl 'tail -3 /root/autodl-tmp/SMVT_AutoDL_Package/trajectories/BIOTIN/BIOTIN_100ns.csv'

# Check GPU
ssh autodl nvidia-smi

# Download results
scp -r autodl:/root/autodl-tmp/SMVT_AutoDL_Package/trajectories/ ./results/

# All done marker
ssh autodl 'cat /root/autodl-tmp/SMVT_AutoDL_Package/ALL_DONE.txt'
```
