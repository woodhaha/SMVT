# colab_md_top3 — Debug Changes

## Issues Fixed (v1 → v2)

| # | Bug | Fix |
|---|-----|-----|
| 1 | Google Drive mount (interactive OAuth, fails on CLI) | Files uploaded directly via `colab upload` to `/content/` |
| 2 | `!cp -r {out_dir} /content/drive/...` — shell `!` doesn't expand Python f-strings | Replaced with `subprocess.run(['cp', '-r', ...])` |
| 3 | `openff-toolkit` unavailable on Colab Python 3.12 | Switched to GAFF2 (`gaff-2.11`) via `openmmforcefields.SystemGenerator` |
| 4 | `MDAnalysis` analysis cell used wrong API (`R.results.time` → correct: `R.results['time']`) | Fixed MDAnalysis API calls |
| 5 | Missing `matplotlib` in pip install | Added to setup cell |

## Colab CLI SSL Fix

Required env vars for colab CLI through WSL with Palantir proxy:
```bash
export SSL_CERT_FILE=/etc/ssl/certs/ca-certificates.crt
export REQUESTS_CA_BUNDLE=/etc/ssl/certs/ca-certificates.crt
```

## Force Field Change

| Original | Debugged |
|----------|----------|
| OpenFF Sage 2.1 (`openff-2.1.0`) | GAFF2 (`gaff-2.11`) |
| Requires `openff-toolkit` | Bundled in `openmmforcefields` |
| `SystemGenerator(..., molecules=[off_mol])` | `SystemGenerator(..., small_molecule_forcefield='gaff-2.11')` |

## Files Needed on Colab VM
- `AF-Q9Y289-F1.pdb` → `/content/`
- `NAFTAZONE_docked.pdbqt` → `/content/`
- `PHENOBARBITAL_docked.pdbqt` → `/content/`
- `ESKETAMINE_docked.pdbqt` → `/content/`

## Run Commands
```bash
# Upload
colab --auth=adc upload -s seer-gastric local_file.ipynb /content/local_file.ipynb

# Execute (always include SSL_CERT_FILE + REQUESTS_CA_BUNDLE)
colab --auth=adc exec -s seer-gastric -f colab_md_top3.ipynb --timeout 86400
```
