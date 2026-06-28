# SMVT MD — Experimental Plan

> Created: 2026-06-28 · Status: Debugging Esketamine pilot

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
