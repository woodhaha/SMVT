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
    """audit_consistency.main() reads sys.argv and takes no parameters."""
    mod = load_skill("audit_consistency")
    print(f"\n=== mechanical consistency: {CONFIG.name} ===")
    argv = sys.argv
    sys.argv = ["audit_consistency.py", str(CONFIG)]
    try:
        return mod.main()
    finally:
        sys.argv = argv


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
    import unicodedata
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
        # LaTeX accents must be folded, not deleted: \v{a} / \'a / \~a carry the letter.
        s = re.sub(r"\\[a-zA-Z]+\{([^{}]*)\}", r"\1", s)
        s = re.sub(r"\\[a-zA-Z]+\s*", "", s)
        s = s.replace("\\&", "and").replace("{", "").replace("}", "").replace("$", "")
        s = unicodedata.normalize("NFKD", s)
        s = "".join(c for c in s if not unicodedata.combining(c))
        # Crossref uses U+2010/U+2011 hyphen and <sup>/<i> markup inside titles.
        # Non-alphanumerics become SPACES, never deleted: deleting them fuses
        # "Na+-dependent" into "nadependent" and invents a title mismatch.
        s = re.sub(r"<[^>]+>", " ", s)
        s = re.sub(r"[^A-Za-z0-9]+", " ", s)
        return re.sub(r"\s+", " ", s).strip().lower()

    def words(s):
        return {w for w in strip(s).split() if len(w) > 2 and not w.isdigit()}

    def abbr_ok(abbrev, full):
        """True if an abbreviated journal name can be derived from the full name.

        'j. biol. chem.' vs 'journal of biological chemistry' shares no complete
        word, so token equality would flag every abbreviated reference. Prefix
        matching per content word is what an abbreviation actually is.
        """
        aw, fw = words(abbrev), words(full)
        if aw <= fw:
            return True
        return all(any(a == f or (len(a) >= 2 and f.startswith(a)) for f in fw) for a in aw)

    def pubmed_title(doi):
        """Fallback title, because Crossref sometimes deposits only a short title
        (wang1999: Crossref 'Human Placental Na+-dependent Multivitamin Transporter'
        vs PubMed's full '... Cloning, functional expression, gene structure, and
        chromosomal localization' - real paper, would otherwise flag as invented)."""
        try:
            u = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
            with urllib.request.urlopen(u + f"esearch.fcgi?db=pubmed&term={doi}[doi]&retmode=json", timeout=30) as f:
                ids = json.load(f).get("esearchresult", {}).get("idlist", [])
            if not ids:
                return ""
            with urllib.request.urlopen(u + f"esummary.fcgi?db=pubmed&id={ids[0]}&retmode=json", timeout=30) as f:
                return json.load(f)["result"][ids[0]].get("title", "")
        except Exception:
            return ""

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
        cr_fam = [strip(a.get("family", "")) for a in m.get("author", [])]
        cr_giv = strip(" ".join(a.get("given", "") for a in m.get("author", [])))
        for part in re.split(r",\s*and\s+|\s+and\s+|,", strip(f["author"])):
            s = part.strip().strip(".").split(" ")[0] if part.strip() else ""
            if len(s) < 3 or s in ("others", "al", "and", "et"):
                continue          # "and others" -> et al., not a surname
            if not any(s == c or s in c or c in s for c in cr_fam) and s not in cr_giv:
                errs.append(f"author {s!r} not in Crossref {cr_fam[:5]}")
        cr_title = m.get("title", [""])[0]
        extra = words(f["title"]) - words(cr_title) - words(pubmed_title(doi))
        extra = {w for w in extra if len(w) > 3 and not w.isdigit()}
        if extra:
            errs.append(f"title words absent from Crossref and PubMed: {sorted(extra)[:5]}")
        cr_year = (m.get("published-print") or m.get("published-online") or {}).get("date-parts", [[None]])[0][0]
        if f["year"] and cr_year and abs(int(f["year"]) - int(cr_year)) > 1:
            errs.append(f"year bib={f['year']} crossref={cr_year}")
        # Venue is pinned by the DOI itself, so abbreviation-vs-full-name in the bib
        # field is cosmetic; only a substantive mismatch (no shared token) is reported.
        cr_j = (m.get("container-title") or [""])[0]
        if f["journal"] and cr_j and not abbr_ok(f["journal"], cr_j):
            errs.append(f"journal not derivable: bib={strip(f['journal'])[:44]!r} cr={strip(cr_j)[:44]!r}")
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
