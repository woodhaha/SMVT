"""AutoDock Vina docking for all SMVT-8 compounds against the SMVT model."""
import os, subprocess, json

ANALYSIS = os.path.join(os.path.dirname(os.path.abspath(__file__)))
LIGANDS = os.path.join(ANALYSIS, 'ligands')
RECEPTOR = os.path.join(ANALYSIS, 'receptor_clean.pdbqt')
VINA_BIN = r'C:\Users\woodh\bin\vina.exe'

# Binding pocket centroid from MM-GBSA minimized complex (BIOTIN LIG centroid)
CENTER = (3.9, 2.0, -7.6)
BOX_SIZE = 22  # Å

COMPOUNDS = {
    'BIOTIN': 'BIOTIN.pdbqt',
    'ESKETAMINE': 'ESKETAMINE.pdbqt',
    'FUROSEMIDE': 'FUROSEMIDE.pdbqt',
    'GABAPENTIN_ENACARBIL': 'GABAPENTIN_ENACARBIL.pdbqt',
    'HYDROMORPHONE': 'HYDROMORPHONE.pdbqt',
    'NAFTAZONE': 'NAFTAZONE.pdbqt',
    'PHENOBARBITAL': 'PHENOBARBITAL.pdbqt',
    'RIBOFLAVIN': 'RIBOFLAVIN.pdbqt',
}

results = {}
for name, lig_pdbqt in COMPOUNDS.items():
    lig_path = os.path.join(LIGANDS, lig_pdbqt)
    out_path = os.path.join(ANALYSIS, f'{name}_vina_out.pdbqt')
    log_path = os.path.join(ANALYSIS, f'{name}_vina_log.txt')

    if not os.path.exists(lig_path):
        print(f'{name}: skip (no ligand pdbqt)')
        continue

    cmd = [
        VINA_BIN,
        '--receptor', RECEPTOR,
        '--ligand', lig_path,
        '--out', out_path,
        '--center_x', str(CENTER[0]),
        '--center_y', str(CENTER[1]),
        '--center_z', str(CENTER[2]),
        '--size_x', str(BOX_SIZE),
        '--size_y', str(BOX_SIZE),
        '--size_z', str(BOX_SIZE),
        '--exhaustiveness', '12',
        '--num_modes', '9',
    ]
    print(f'{name}: docking...', end=' ', flush=True)
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    print(f'done (rc={r.returncode})', flush=True)
    # Save log from stderr
    with open(log_path, 'w') as f:
        f.write(r.stderr)
        f.write(r.stdout)

    scores = []
    if r.returncode == 0:
        for ln in r.stdout.split('\n'):
            parts = ln.strip().split()
            if len(parts) >= 2 and parts[0] == '1':
                try:
                    scores.append(float(parts[1]))
                except: pass
    results[name] = {
        'best_score': scores[0] if scores else None,
    }
    if scores:
        print(f'  best: {scores[0]:.3f} kcal/mol')

# Summary
print(f'\n{"="*60}')
print(f'  VINA DOCKING SUMMARY')
print(f'{"="*60}')
print(f'  {"Compound":30s}  {"Score (kcal/mol)":>15s}')
print(f'  {"-"*48}')
for name, r in sorted(results.items(), key=lambda x: x[1].get('best_score', 0) if x[1].get('best_score') else 0):
    sc = r.get('best_score')
    sc_str = f'{sc:+.1f}' if sc else 'FAILED'
    print(f'  {name:30s}  {sc_str:>15s}')

json.dump(results, open(os.path.join(ANALYSIS, 'vina_results.json'), 'w'), indent=2)
print(f'\nSaved to vina_results.json')
