#!/usr/bin/env python3
import re, glob
tex = open(r'D:\Researching\SMVT\04_Manuscript\SMVT_manuscript.tex', encoding='utf-8').read()
bib = open(r'D:\Researching\SMVT\04_Manuscript\references.bib', encoding='utf-8').read()
bbl = open(r'D:\Researching\SMVT\04_Manuscript\SMVT_manuscript.bbl', encoding='utf-8').read()
cite = set()
for m in re.finditer(r'\\cite\{([^}]*)\}', tex):
    for k in m.group(1).split(','):
        cite.add(k.strip())
bibkeys = set(re.findall(r'@\w+\{([^,]+),', bib))
bblkeys = set(re.findall(r'\\bibitem\{([^}]+)\}', bbl))
print('cite count:', len(cite), '| bib keys:', len(bibkeys), '| bbl keys:', len(bblkeys))
print('cited but MISSING from bib:', sorted(cite - bibkeys))
print('bib entries NEVER cited:', sorted(bibkeys - cite))
print('bib vs bbl mismatch:', sorted(bibkeys ^ bblkeys))
# abstract word count
m = re.search(r'\\begin\{abstract\}(.*?)\\end\{abstract\}', tex, re.S)
body = m.group(1)
body = re.sub(r'\\[a-zA-Z]+\*?', ' ', body)
body = body.replace('{','').replace('}','').replace('$','').replace('~',' ')
words = re.findall(r'[A-Za-z0-9%.\-/]+', body)
print('abstract words:', len(words))
