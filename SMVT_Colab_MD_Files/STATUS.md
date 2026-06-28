# SMVT MD on Colab — Debug Status

## What Works
- ✅ Colab GPU T4 session (`seer-gastric`) running
- ✅ OpenMM 8.5.2 on OpenCL (GPU)
- ✅ All files uploaded to `/content/`
- ✅ Protein preparation (PDBFixer: add missing atoms, H at pH 7.4)
- ✅ Vina pose extraction (19, 19, 17 atoms)
- ✅ RDKit SMILES → 3D PDB generation
- ✅ Colab CLI with SSL fix (`SSL_CERT_FILE` + `REQUESTS_CA_BUNDLE`)

## What Doesn't Work
- ❌ `openff-toolkit` NOT available on Colab Python 3.12 pip
- ❌ `openff-toolkit` is REQUIRED by both GAFF and SMIRNOFF ligand parameterization
- ❌ `openmmforcefields.SystemGenerator` cannot work without it
- ❌ Manual residue template generation is complex (GAFF atom typing needed)
- ❌ Esketamine has Cl atoms — needs special GAFF parameters

## The Blocker
Colab's pip index doesn't serve `openff-toolkit` for Python 3.12.
Without it, OpenMM can't parameterize non-standard small molecule residues.
Both `openmmforcefields` GAFF and SMIRNOFF backends require `import openff.toolkit`.

## Solution: Local Parameterization + Upload

1. **Generate ligand parameters locally** (Windows, with conda openff-toolkit):
```bash
conda install -c conda-forge openff-toolkit
python generate_ligand_params.py
```

2. **Upload parameter files to colab**:
```bash
colab --auth=adc upload -s seer-gastric ligand_params.xml /content/
```

3. **Run MD with ForceField + custom templates**.

## Files Saved
- `colab_md_top3_debugged.ipynb` — debugged notebook (GAFF2, no Drive, f-string fixes)
- `README_debug.md` — debug change log
- `STATUS.md` — this file

## Colab CLI SSL Fix (Critical)
```bash
export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
export REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
```
Both are REQUIRED for colab CLI through WSL with Palantir proxy.
