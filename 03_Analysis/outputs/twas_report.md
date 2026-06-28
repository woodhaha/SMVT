# SMVT TWAS — SLC5A6 → Colorectal Cancer Risk

**Date**: 2026-06-23
**Method**: S-PrediXcan / TWAS framework

## Data Sources

- **eQTL**: GTEx v8 (colon sigmoid n=318, colon transverse n=368, whole blood n=670)
- **GWAS**: FinnGen CRC / Huyghe et al. 2019 (n~120K)

## Results

### SLC5A6 eQTL Status

**No significant single-tissue cis-eQTLs were found for SLC5A6 in any GTEx v8 colon tissue.**


### TWAS Power Estimation

- cis-h² (estimated): 0.10
- N_eff: 317
- Expected Z_TWAS: 0.398
- TWAS power (Bonferroni α=2.5e-6): 0.0%

### Colocalization

- Nearest CRC GWAS locus: rs35360328 (2p24.1), ~2.2 Mb upstream of SLC5A6
- No evidence for shared causal variant (H4 < 0.1)

## Conclusion

**TWAS is underpowered** for SLC5A6 → CRC at current sample sizes. SLC5A6 lacks strong cis-eQTL instruments in colon tissue, and its expression heritability in colon is likely modest. The gene is expression-driven rather than genetically regulated, which is consistent with our finding that SMVT is an expression-driven (not mutation-driven) target.

### Alternative strategies recommended:
1. **Larger eQTL reference**: Wait for GTEx v9/v10 with larger colon sample sizes
2. **Multi-tissue TWAS**: Use UTMOST or MetaXcan to borrow strength across tissues
3. **Focus on functional validation**: The docking/KO/survival evidence is stronger than genetic evidence for this target
