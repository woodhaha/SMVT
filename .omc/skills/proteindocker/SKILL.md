---
name: proteindocker
description: Use when designing peptide therapeutics de novo, when a target has no typical small-molecule pocket (PPI interfaces, flat surfaces), when virtual screening hits are exhausted, or when the user asks about peptide design, de novo peptide generation, or ProteinDocker workflows.
---

# ProteinDocker — Peptide De Novo Design

## Overview

Virtual screening picks from existing libraries. ProteinDocker **creates peptides from scratch** — generative algorithms design novel sequences and conformations tailored to a target's 3D structure, unconstrained by any molecular database.

## When to Use

- Target has **no druggable small-molecule pocket** (PPI interfaces, flat surfaces)
- Virtual screening exhausted without satisfactory hits
- Need high-specificity binders with large contact interfaces
- Designing peptide inhibitors, probes, or functional peptides

**NOT for:** routine small-molecule docking (use Vina/Glide), antibody design, or cases where an existing library suffices.

## Paradigm

| Virtual Screening | ProteinDocker (De Novo) |
|---|---|
| Screens existing molecules | Generates novel sequences |
| Diversity capped by library | Library-free — infinite space |
| "Find the key in the warehouse" | "Make the key for the lock" |

## 5-Stage Pipeline

### 1. Target Preprocessing
- Parse target protein 3D structure
- Identify and isolate the binding interface
- Clean (remove waters/heteroatoms), repair missing loops, energy-minimize

### 2. De Novo Peptide Generation
- Generative algorithm explores peptide length, amino acid composition, and 3D conformation
- Output: novel sequences with target-specific binding poses
- No database constraint — every peptide is new

### 3. Structure Refinement
- Fix unreasonable conformations and atomic clashes
- Optimize intermolecular interactions (H-bonds, salt bridges, hydrophobic packing)
- Improve peptide-target complementarity

### 4. Multi-Dimensional Screening
Rank candidates on:
- **Interface completeness** — does the peptide cover the target region?
- **Conformational plausibility** — are geometries physically reasonable?
- **Binding energy** — scoring function ranking
- **System stability** — coarse energetic assessment

### 5. MD Validation
- Molecular dynamics simulation to verify binding stability under thermal motion
- MM/GBSA binding free energy calculation
- Output: a small set of candidates ready for wet-lab synthesis and testing

## Peptide Properties

### Advantages
- **High specificity** — multi-residue contacts anchor to broad, flat interfaces; low off-target risk
- **Highly engineerable** — tune length, substitute residues, cyclize, add disulfide bonds or non-natural amino acids
- **Cost-effective vs antibodies** — smaller, simpler, mature chemical synthesis

### Limitations (computational predictions ≠ experimental reality)
- Linear peptides degraded by proteases (short plasma half-life)
- Poor membrane penetration
- Potential solubility/aggregation/toxicity issues
- **ALL candidates require experimental validation** — computation narrows the search space, it doesn't replace synthesis and assay

## Common Mistakes

- **Skipping target prep** — dirty structures (clashes, missing residues) poison generation
- **Over-relying on scores** — a top-ranked peptide can still fail in MD; always run dynamics
- **Ignoring developability** — a perfect binder that aggregates or degrades instantly is not a drug
- **Treating computation as the endpoint** — wet-lab validation is non-negotiable

## References

Source: 从"库里找分子"到"按需造分子"：多肽药物设计迎来全新范式 (WeChat, mp.weixin.qq.com)
