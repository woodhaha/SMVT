"""ADMET prediction via RDKit — run in smvt-md conda env"""
import json
import sys
import os

sys.path.insert(0, r'C:\anaconda3\envs\smvt-md\Lib\site-packages')

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED

ANALYSIS = r'D:\Researching\SMVT\MM-GBSA Analysis'

with open(os.path.join(ANALYSIS, 'compounds.json')) as f:
    data = json.load(f)

print(f'{"Compound":20s} {"MW":>7s} {"LogP":>6s} {"HBD":>4s} {"HBA":>4s} {"RB":>4s} {"TPSA":>6s} {"Ro5":>4s} {"QED":>5s}')
print('-'*65)

results = {}
for c in data['compounds']:
    mol = Chem.MolFromSmiles(c['smiles'])
    if mol is None:
        print(f'{c["id"]:20s}  INVALID SMILES')
        continue
    mw = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rb = Lipinski.NumRotatableBonds(mol)
    tpsa = Descriptors.TPSA(mol)
    ro5 = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    qed = QED.qed(mol)

    # Simple BBB heuristic
    if logp > 1.5 and logp < 3.5 and tpsa < 90 and mw < 450:
        bbb = 'Yes'
    elif logp > 0.5 and logp < 4.5 and tpsa < 120 and mw < 500:
        bbb = 'Borderline'
    else:
        bbb = 'No'

    print(f'{c["id"]:20s} {mw:6.1f} {logp:5.2f} {hbd:3d} {hba:3d} {rb:3d} {tpsa:5.1f} {ro5:3d} {qed:.3f}')

    results[c['id']] = {
        'MW': round(mw, 1), 'LogP': round(logp, 2),
        'HBD': hbd, 'HBA': hba, 'RotBonds': rb,
        'TPSA': round(tpsa, 1), 'Ro5_violations': ro5,
        'QED': round(qed, 3), 'BBB_heuristic': bbb,
    }

with open(os.path.join(ANALYSIS, 'admet_predictions.json'), 'w') as f:
    json.dump(results, f, indent=2)

print(f'\n=== ADMET SUMMARY ===')
print(f'{"Compound":20s} {"MW":>6s} {"LogP":>6s} {"TPSA":>5s} {"BBB":>8s} {"Ro5":>4s} {"QED":>5s}')
print('-'*55)
for rid, r in results.items():
    print(f'{rid:20s} {r["MW"]:5.1f} {r["LogP"]:5.2f} {r["TPSA"]:5.0f} {r["BBB_heuristic"]:>8s} {r["Ro5_violations"]:3d} {r["QED"]:.3f}')

print(f'\nSaved: admet_predictions.json')
