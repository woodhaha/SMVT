# SMVT — Master Plan & Results

> Updated: 2026-06-28 · Status: Esketamine MD debugging, screening complete

---

## Virtual Screening Results (completed)

### Pipeline: ZINC 230M → FDA/ChEMBL 3,300 → ML (AUC=0.888) → 500 docked → 8 Elite Hits

### Top 8 Elite Hits (ΔG < −8.0 kcal/mol)

| Rank | Compound | ΔG | Class | Novel? |
|:--:|----------|:--:|-------|:---:|
| 1 | Hydromorphone | −8.58 | Opioid | ✅ First opioid SMVT ligand |
| 2 | Furosemide | −8.36 | Sulfonamide | ✅ New scaffold |
| 3 | Naftazone | −8.34 | Naphthoquinone | — |
| 4 | Phenobarbital | −8.30 | Barbiturate | 100% class hit rate |
| 5 | Pentobarbital | −8.18 | Barbiturate | |
| 6 | Diclofenac | −8.07 | NSAID | NSAID-SMVT axis |
| 7 | Carprofen | −8.04 | NSAID | |
| 8 | Butabarbital | −8.02 | Barbiturate | |

### Chemical Family Hit Rates

| Family | Hit/Tested | Rate | Key Insight |
|--------|:---:|:---:|------|
| Barbiturate | 8/8 | **100%** | Barbituric acid = biotin ureido mimic |
| NSAID | 6/12 | 50% | Confirms NSAID-SMVT transport axis |
| Opioid | 2/5 | 40% | Morphinan = new SMVT ligand class |
| Sulfonamide | 1/6 | 17% | Furosemide hit |
| Vitamin/Cofactor | 2/8 | 25% | Biotin baseline −6.76 |

### Pharmacophore SAR Rules
1. Cyclic ureide/carboxamide core = key pharmacophore
2. Carboxyl NOT essential (barbiturates lack it → expands chemical space)
3. Aromatic ring enhances affinity (hydrophobic contact)
4. ≥2 H-bond acceptors (C=O) for optimal binding
5. Halogen (Cl) tolerated, adds specificity

### Controls
- Biotin (−6.76): natural substrate reference
- Gabapentin enacarbil (−6.63): FDA-approved SMVT prodrug, positive control
- Riboflavin (−0.01): non-binding vitamin, negative control

---

## Experiment Design

**Target**: SMVT (SLC5A6) · Receptor: AF-Q9Y289-F1 (AlphaFold, pLDDT=79.4)
**Protocol**: Conflict check → GAFF template → Staged min → NVT→NPT → 100ns production
**Analysis**: RMSD / RMSF / H-bond / MM-GBSA / SASA

## Compound Matrix (7 + 1 optional)

### Test Compounds (Top 4 from virtual screening)

| # | Compound | ΔG | Class | Status |
|:--:|----------|:--:|-------|--------|
| 1 | Hydromorphone | −8.58 | Opioid | ⬜ |
| 2 | Furosemide | −8.36 | Sulfonamide (Cl+S) | ⬜ |
| 3 | Naftazone | −8.34 | Naphthoquinone | ⬜ |
| 4 | Phenobarbital | −8.30 | Barbiturate | ⬜ |

### Controls

| # | Compound | ΔG | Class | Role | Status |
|:--:|----------|:--:|-------|------|--------|
| 5 | Biotin | −6.76 | Vitamin B7 | 🔵 Natural substrate | ⬜ |
| 6 | Gabapentin enacarbil | −6.63 | Prodrug | 🟢 Positive (FDA) | ⬜ |
| 7 | Riboflavin | −0.01 | Vitamin B2 | 🔴 Negative (no binding) | ⬜ |

### Pilot (debugging probe)

| # | Compound | ΔG | Class | Status |
|:--:|----------|:--:|-------|--------|
| 8* | Esketamine | −7.58 | Arylcyclohexylamine (Cl) | 🔄 Staged min [9b] ✓ |

## Per-Compound Workflow

```
1. analyze_clashes.py        → Check Vina pose for steric clashes
2. gen_gaff_v3.py            → Generate GAFF 2.11 template (RDKit + MMFF94)
3. md_v7_staged_min.py       → Staged min → NVT 50→300K → NPT → Production
4. analyze_md.py             → RMSD/RMSF/H-bond/MM-GBSA/SASA
```

## Root Cause (from Esketamine pilot)

Vina docking pose has **steric clashes** with protein (shortest 2.08Å = 61% vdW).
→ NVT crashes at 200K as kinetic energy overcomes vdW barrier.
→ **Fix**: 4-stage minimization (all restrained → free protein → free ligand → all free)
→ Verified: [9b] PE dropped from 7.6T to −885K.

## Key Files

| File | Purpose |
|------|---------|
| `SMVT_Colab_MD_Files/gen_gaff_v3.py` | GAFF template generator (RDKit + MMFF94) |
| `SMVT_Colab_MD_Files/gen_gaff_templates.py` | GAFF template generator (Gasteiger fallback) |
| `SMVT_Colab_MD_Files/analyze_clashes.py` | Vina pose clash detector |
| `SMVT_Colab_MD_Files/md_v7_staged_min.py` | Staged MD protocol |
| `SMVT_Colab_MD_Files/analyze_md.py` | Post-MD analysis |
| `SMVT_MD_Package/README.md` | Original experiment spec |
| `SMVT_Colab_MD_Files/ligand_params_v3/` | MMFF94 templates (production) |
| `SMVT_Colab_MD_Files/ligand_params/` | Gasteiger templates (fallback) |

## Execution Strategy: Colab Free Tier (Zero Cost)

Rationale: 阿里云竞价实例 ¥1,000-2,000 vs Colab 免费。省下的钱够 Colab Pro+ 一年。

| Option | 7×100ns | Cost | Speed |
|--------|:---:|------|------|
| Colab T4 free | ~30 days (serial) | ¥0 | Slow |
| Colab Pro+ | ~14 days | $50/mo | Medium |
| Ali V100 spot | ~7 days | ¥1,000-1,500 | Fast |
| Ali T4 spot | ~21 days | ¥1,500-2,000 | Slow |

**Plan**: Colab free, 1 compound/day. Add `CheckpointReporter(10ns)` so interrupted runs resume.

### Daily Run

```bash
colab new --gpu T4 -s smvt-md
colab upload AF-Q9Y289-F1.pdb + {COMPOUND}_docked.pdbqt + ligand_params/{COMPOUND}_template.xml
colab install openmm pdbfixer rdkit openmmforcefields
colab exec -f md_v7_staged_min.py --timeout 43200  # 12h timeout
# Download results → analyze → next compound
```

### Run Order (1/day)

1. Esketamine (pilot, finish NVT→production) — verify protocol
2. Biotin (reference) — validate binding site
3. Phenobarbital — simplest test hit
4. Naftazone
5. Furosemide (has Cl+S)
6. Hydromorphone (best hit)
7. Gabapentin enacarbil (positive control)
8. Riboflavin (negative control)

## Resume Command

```bash
colab new --gpu T4
colab upload AF-Q9Y289-F1.pdb + *_docked.pdbqt + ligand_params/*.xml
colab install openmm pdbfixer rdkit openmmforcefields
colab exec -f md_v7_staged_min.py --timeout 7200
```

## Expected Outcomes

- **Hits (#1-4)**: RMSD < 3Å → stable binding confirmation
- **Biotin (#5)**: RMSD 2-4Å → natural substrate dynamics
- **Gabapentin (#6)**: Higher RMSD → transient transporter interaction
- **Riboflavin (#7)**: Dissociation → validates screening specificity
