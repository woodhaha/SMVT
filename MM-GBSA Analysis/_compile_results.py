"""Compile MM-GBSA results summary."""
import json, os

analysis = os.path.join(os.path.dirname(os.path.abspath(__file__)))
comps = ['NAFTAZONE','BIOTIN','ESKETAMINE','FUROSEMIDE',
         'GABAPENTIN_ENACARBIL','HYDROMORPHONE','PHENOBARBITAL','RIBOFLAVIN']

results = {}
for c in comps:
    fn = os.path.join(analysis, c, 'mmgbsa_results.json')
    if os.path.exists(fn):
        with open(fn) as f:
            results[c] = json.load(f)

print('=' * 72)
print('  MM-GBSA BINDING FREE ENERGY -- SMVT-8 COMPOUNDS')
print('  OpenMM amber14 + GB(OBC2) + GAFF2.11')
print('  Single-structure (last frame 100ns MD, restrained min)')
print('=' * 72)
print(f'  {"Compound":28s}  {"dG(kJ/mol)":>10s}  {"dG(kcal/mol)":>10s}  {"Lig":>4s}  {"Rank":>4s}')
print(f'  {"-"*28}  {"-"*10}  {"-"*10}  {"-"*4}  {"-"*4}')

# Sort by dG_kcal (most negative = best binder)
sorted_compounds = sorted(results.items(), key=lambda x: x[1]['dG_kcal'])
flags = {
    'FUROSEMIDE': ' INVAL',
    'RIBOFLAVIN': ' NEG',
    'BIOTIN': ' REF',
    'GABAPENTIN_ENACARBIL': ' PRODRUG',
}

for rank, (c, r) in enumerate(sorted_compounds, 1):
    flag = flags.get(c, '')
    print(f'  {c:28s}  {r["dG_kJ"]:8.1f} kJ  {r["dG_kcal"]:8.2f} kcal  {r["n_lig_atoms"]:3d}  #{rank}{flag}')

print()
print('  Rankings (by dG):')
for rank, (c, r) in enumerate(sorted_compounds, 1):
    flag = flags.get(c, '')
    flag_txt = f' [{flag.strip()}]' if flag else ''
    print(f'  #{rank}: {c}{flag_txt} ({r["dG_kcal"]:+.2f} kcal/mol)')

print()
print('  Notes:')
print('  - FUROSEMIDE: template is a fragment (missing furan ring). dG invalid.')
print('  - RIBOFLAVIN: negative control, unexpected strong binding predicted.')
print('  - Single-structure MM-GBSA, no entropy (-TdS) term included.')
print(f'  - Results: {analysis}/<COMPOUND>/mmgbsa_results.json')
