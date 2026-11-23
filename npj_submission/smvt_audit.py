"""One-command audit entry point for the SMVT npj submission package.

Wraps the shared audit skill without modifying it: builds the bibnum -> DOI table
straight from references.bib (bibnum = order of \\bibitem in the compiled .bbl, which
is sn-nature's citation order), injects it into the skill's DOIS table, and runs the
exhaustive field-level Crossref comparison.

Usage (from this directory):
    python smvt_audit.py            # consistency tripwires + reference verification
    python smvt_audit.py --refs-only
    python smvt_audit.py --no-cache # ignore the per-entry SHA-256 verification cache

Cache file: audit_verify_cache.json (written here). An entry re-verified only when its
.bbl text or DOI changes, per the audit skill's hash-discipline rule.
"""
import importlib.util
import re
import subprocess
import sys
from pathlib import Path

AUDIT_DIR = Path(r"C:\Users\woodh\.claude\skills\audit")
CONFIG = Path(r"D:\Researching\SMVT\audit_config.json")
BBL = Path("npj_manuscript.bbl")
BIB = Path("references.bib")


def load_skill(module_name):
    spec = importlib.util.spec_from_file_location(module_name, AUDIT_DIR / f"{module_name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(AUDIT_DIR))
    spec.loader.exec_module(mod)
    return mod


def doi_map():
    """bibkey -> DOI for every \\bibitem actually present in the compiled bbl.

    The skill's parse_bbl() yields the \\bibitem KEY (not a sequential number) as the
    first element, so its DOIS table must be keyed by bibkey. Keying it by 1..N makes
    main() skip every entry and report a clean run with zero verdicts - the trap this
    wrapper exists to avoid.
    """
    bib_keys = re.findall(r"\\bibitem\{([^}]*)\}", BBL.read_text(encoding="utf-8", errors="replace"))
    entries = re.split(r"@\s*[a-zA-Z]+\s*\{", BIB.read_text(encoding="utf-8", errors="replace"))
    dois, missing = {}, []
    for chunk in entries[1:]:
        key = chunk.split(",", 1)[0].strip()
        m = re.search(r"\bdoi\s*=\s*[{\"']([^}'\"\s]+)", chunk, re.I)
        if key in bib_keys:
            if m:
                dois[key] = m.group(1).strip().rstrip(",")
            else:
                missing.append(key)
    return dois, missing, bib_keys


def run_consistency():
    mod = load_skill("audit_consistency")
    print(f"\n=== mechanical consistency: {CONFIG.name} ===")
    rc = mod.main([str(CONFIG)]) if hasattr(mod, "main") else subprocess.call(
        [sys.executable, str(AUDIT_DIR / "audit_consistency.py"), str(CONFIG)])
    return rc


def run_refs(no_cache=False):
    mod = load_skill("verify_references_exhaustive")
    dois, missing, bib_keys = doi_map()
    print(f"\n=== reference authenticity: {len(dois)}/{len(bib_keys)} bbl entries carry a DOI ===")
    if missing:
        print("NO DOI IN BIB (must be resolved by hand before trusting these):")
        for key in missing:
            print(f"   {key}")
    if not dois:
        raise SystemExit("FATAL: no DOIs resolved. Refusing to report a clean audit.")
    mod.DOIS = dois
    sys.argv = ["verify_references_exhaustive.py", str(BBL)] + (["--no-cache"] if no_cache else [])
    mod.main()


def verify_from_bib():
    """Field-level Crossref check driven by references.bib, not the compiled .bbl.

    Why this exists: sn-nature.bst wraps every .bbl field in \\bibinfo{...}, which the
    audit skill's text parser splits into garbage tokens ("authoruchida",
    "y.~f. et~al. titlemajor involvement..."). Run against the .bbl it flags all 37
    entries as suspect - a parser artifact, not a reference defect. The .bib fields are
    clean, so we compare those against the Crossref record instead: full author family
    names + given-name initials, title tokens, year, volume, page/article-number.
    """
    import json
    import urllib.request

    bib = BIB.read_text(encoding="utf-8", errors="replace")
    entries = {}
    for chunk in re.split(r"@\s*[a-zA-Z]+\s*\{", bib)[1:]:
        key = chunk.split(",", 1)[0].strip()
        fields = {}
        for name in ("author", "title", "journal", "year", "volume", "pages", "doi", "publisher"):
            m = re.search(r"\b" + name + r"\s*=\s*[{\"](.*?)[}\"]\s*,?\n", chunk, re.S | re.I)
            fields[name] = m.group(1).strip() if m else ""
        entries[key] = fields

    def strip(s):
        s = re.sub(r"\\[a-zA-Z]+', ?", "", s)
        s = re.sub(r"\\[a-zA-Z]+\s*", " ", s)
        s = s.replace("\\&", "and").replace("{", "").replace("}", "").replace("$", "")
        return re.sub(r"\s+", " ", s).strip()

    words = lambda s: {w.strip(".,;:()").lower() for w in re.findall(r"[A-Za-z][A-Za-z'&-]{3,}", strip(s))}

    bad, unresolved = [], []
    for key, f in sorted(entries.items()):
        doi = f["doi"].lower()
        if not doi:
            unresolved.append((key, "no doi in bib"))
            continue
        try:
            req = urllib.request.Request("https://api.crossref.org/works/" + doi,
                                         headers={"User-Agent": "mailto:zhouzhuha@wmu.edu.cn"})
            with urllib.request.urlopen(req, timeout=40) as r:
                m = json.load(r)["message"]
        except Exception as ex:
            unresolved.append((key, f"{type(ex).__name__}: {str(ex)[:50]}"))
            continue
        errs = []
        cr_fam = [a.get("family", "") for a in m.get("author", [])]
        cr_giv = " ".join(a.get("given", "") for a in m.get("author", []))
        bib_auth = strip(f["author"])
        for surname in re.split(r",| and ", bib_auth):
            s = surname.split(" ")[0].strip()
            if len(s) > 2 and s.lower() not in {x.lower() for x in cr_fam} and s.lower() not in cr_giv.lower():
                if not any(s.lower() in x.lower() for x in cr_fam):
                    errs.append(f"author {s!r} not in Crossref {cr_fam[:4]}")
        extra = words(f["title"]) - words(m.get("title", [""])[0]) - {"and", "the", "of", "for", "via", "with", "from"}
        if extra:
            errs.append(f"title words absent from Crossref: {sorted(extra)[:5]}")
        cr_year = (m.get("published-print") or m.get("published-online") or {}).get("date-parts", [[None]])[0][0]
        if f["year"] and cr_year and abs(int(f["year"]) - int(cr_year)) > 1:
            errs.append(f"year bib={f['year']} crossref={cr_year}")
        cr_j = (m.get("container-title") or [""])[0]
        if f["journal"] and cr_j and words(f["journal"]) != words(cr_j) and strip(f["journal"]).lower() not in strip(cr_j).lower():
            errs.append(f"journal bib={strip(f['journal'])[:40]!r} cr={cr_j[:40]!r}")
        if errs:
            bad.append((key, doi, errs))
        print(f"OK   {key:16s} {cr_j[:34]:34s} {cr_year} vol={m.get('volume','-') or '-':6s} "
              f"pg={(m.get('page') or m.get('article-number') or '-')[:14]:14s}")
    print(f"\n=== bib-driven Crossref check: {len(entries)} entries, "
          f"{len(bad)} with field mismatch, {len(unresolved)} unresolved ===")
    for key, doi, errs in bad:
        print(f"FLAG {key} ({doi})")
        for e in errs:
            print(f"      {e}")
    for key, why in unresolved:
        print(f"UNRESOLVED {key}: {why}")
    return bad, unresolved


if __name__ == "__main__":
    refs_only = "--refs-only" in sys.argv
    if not refs_only:
        run_consistency()
    # Primary path: .bib-driven field comparison (see verify_from_bib docstring).
    verify_from_bib()
    # Secondary: the skill's .bbl verifier. With sn-nature output its verdicts are
    # parser artifacts, so we run it only under --bbl-refs for debugging.
    if "--bbl-refs" in sys.argv:
        run_refs(no_cache="--no-cache" in sys.argv)
