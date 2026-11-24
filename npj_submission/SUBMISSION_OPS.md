# SNAPP submission operations — npj Systems Biology and Applications

Run `python stage_snapp.py` from this directory to (re)build the upload set. It reads
the live files here, so re-run it on the day you submit — the staged copy will
otherwise carry the previous day's cover-letter date and any figure rendered since.

Staged to `C:\Users\woodh\Documents\_snapp_npj\` (must be under the workspace root:
browser-automation `upload_file` cannot read anything under `D:\Researching`):

| Artifact | What SNAPP wants it for |
|---|---|
| `manuscript_source.zip` | revised/original manuscript — tex + `references.bib` + `npj_manuscript.bbl` + `sn-jnl.cls` + `figures/*.pdf`. The `.bbl` and the class file must be inside: the server compiles LaTeX but does not run BibTeX |
| `main_figures/Fig1..9.tiff` | the figure field — accepts only `.jpeg/.tiff/.eps`, **not PDF**. All nine are 600 dpi LZW |
| `cover_letter.txt` | the cover-letter field (plain text; pandoc-to-PDF trips on LaTeX math, and markdown `#`/`**` markers would reach the editor) |
| `MANIFEST.txt` | sizes + MD5 of every staged file — verify the upload against this instead of eyeballing the review page |

## Upload order and known traps

1. **Files** → manuscript zip, then the nine TIFFs **serially**. Parallel
   `upload_file` calls against the same button silently drop files (an earlier
   submission showed "9 FILES ADDED" where 12 were sent). Check the count on the
   Review page before continuing.
2. **Details** → title/abstract are auto-extracted from the source; verify both.
   Cover letter has its own upload box (`id=cover_letter`), not Related files.
3. **Authors** → every author needs an e-mail. **`Jiancheng Sun` is still missing his
   e-mail and ORCID** (`authors.md` marks it; do not invent one). The auto-parser
   mangles name strings of the form `Jiancheng Sun, MD, Attending Surgeon` into
   first="Attending"/last="Surgeon" — check each author after import and fix with the
   React native-setter path, not `fill`.
4. **Declarations** → OA + APC agreement, COI is an on-page radio (no signed PDF),
   data availability textarea.
5. **Review** → preprint question: this manuscript has never been submitted anywhere
   (CBM/CSBJ/PLOS CB were all prepared and abandoned before submission), so no Research
   Square preprint should exist for it — confirm on the preprint list before answering,
   then tick terms and submit.

## Before clicking Submit

- `python verify_package.py` → 33/33 (recomputes every number from source data)
- `python smvt_audit.py` → tripwires pass + 37 references Crossref-clean
- Corresponding ORCID must be `0009-0005-1818-4789` (the Springer-bound one); the
  duplicate `0000-0002-4523-2943` merge is still pending — do not enter it.
- APC budget: npj SBA is fully OA, **$3,490 ex-VAT**, priced at acceptance, no China
  waiver. Confirm funding/waiver status with the co-authors beforehand.
