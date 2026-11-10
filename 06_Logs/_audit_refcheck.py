#!/usr/bin/env python3
"""One-off SMVT reference authenticity check (Crossref). Alphanumeric bibitem keys."""
import re, json, sys, time, urllib.request, urllib.parse
from pathlib import Path

UA = {'User-Agent': 'audit-ref-check/1.0 (mailto:zhouzhuha@wmu.edu.cn)'}
BBL = r"D:\Researching\SMVT\04_Manuscript\SMVT_manuscript.bbl"

def strip(s):
    s = re.sub(r'\$.*?\$', ' ', s)               # math
    s = re.sub(r'\\[a-zA-Z]+\s?', ' ', s)        # \em \newblock \" etc
    s = s.replace('{', '').replace('}', '').replace('~', ' ')
    s = re.sub(r'\s+', ' ', s).strip()
    return s

def parse(path):
    txt = Path(path).read_text(encoding='utf-8')
    refs = []
    for m in re.finditer(r'\\bibitem\{([^}]+)\}\s*(.*?)(?=\\bibitem\{|\s*\\end\{thebibliography\})', txt, re.S):
        key, body = m.group(1), m.group(2)
        parts = body.split(r'\newblock')
        authors = strip(parts[0]).rstrip('.')
        title = strip(parts[1]).strip(' .') if len(parts) > 1 else ''
        rest = parts[2] if len(parts) > 2 else ''
        jm = re.search(r'em\{([^}]*)\}', rest)
        journal = strip(jm.group(1)) if jm else ''
        ym = re.search(r'\b(19|20)\d{2}\b', strip(rest))
        year = int(ym.group()) if ym else None
        # first author surname for fallback
        first = authors.split(',')[0].split()
        surname = first[-1] if first else ''
        refs.append(dict(key=key, title=title, journal=journal, year=year, authors=authors, surname=surname))
    return refs

def crossref(**params):
    q = urllib.parse.urlencode(params)
    url = f"https://api.crossref.org/works?rows=3&select=title,container-title,published,DOI,volume,issue,page,article-number,author&{q}"
    req = urllib.request.Request(url, headers=UA)
    with urllib.request.urlopen(req, timeout=25) as r:
        return json.load(r).get('message', {}).get('items', [])

def jname(it):
    return (it.get('container-title') or ['?'])[0]

def run():
    refs = parse(BBL)
    print(f"# {len(refs)} refs\n")
    for i, e in enumerate(refs, 1):
        items = crossref(**{'query.title': e['title']})
        # fallback: add author surname to bibliographic query if no decent match
        if not items and e['surname']:
            items = crossref(**{'query.bibliographic': f"{e['surname']} {e['title']}"})
        if not items:
            print(f"[{i}] {e['key']} *** NO MATCH ***  title={e['title'][:60]}")
            time.sleep(0.2); continue
        top = items[0]
        cy = (top.get('published') or {}).get('date-parts', [[None]])[0][0]
        cj = jname(top)[:45]
        doi = top.get('DOI','')
        ct = (top.get('title') or ['?'])[0][:60]
        ca = '; '.join(f"{a.get('family','')}" for a in (top.get('author') or [])[:4])
        flags = []
        if e['year'] and cy and abs(e['year']-cy) > 1:
            flags.append(f"YEAR bbl={e['year']} cr={cy}")
        if e['journal'] and cj and e['journal'].lower().split(':')[0] not in cj.lower() and cj.lower().split(':')[0] not in e['journal'].lower():
            flags.append(f"JRN bbl='{e['journal']}' cr='{cj}'")
        flag = (' *** ' + ' | '.join(flags) + ' ***') if flags else ''
        print(f"[{i}] {e['key']} {cy} | {cj} | {doi}{flag}")
        print(f"      bbl: {e['journal']} {e['year']} | cr: {ct}")
        if flags:
            print(f"      authors: {ca}")
        time.sleep(0.2)

if __name__ == '__main__':
    run()
