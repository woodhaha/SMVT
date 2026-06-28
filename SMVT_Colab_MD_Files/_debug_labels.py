"""Quick debug: inspect ff.label_molecules output structure."""
from openff.toolkit import Molecule, Topology
from openff.toolkit.typing.engines.smirnoff import ForceField
import json

mol = Molecule.from_smiles('CCC1(C(=O)NC(=O)NC1=O)C2=CC=CC=C2', allow_undefined_stereo=True)
mol.generate_conformers(n_conformers=1)
ff = ForceField('openff-2.1.0.offxml')
labels = ff.label_molecules(mol.to_topology())[0]

print("=== Label keys:", list(labels.keys()))
for k, v in labels.items():
    if isinstance(v, list):
        print(f"\n{k}: list[{len(v)}]")
        if len(v) > 0:
            item = v[0]
            print(f"  type: {type(item).__name__}")
            if isinstance(item, dict):
                for k2, v2 in list(item.items())[:5]:
                    print(f"  {k2}: {v2}")
            elif isinstance(item, (list, tuple)):
                print(f"  content: {item[:3]}")
    elif isinstance(v, dict):
        print(f"\n{k}: dict[{len(v)}]")
        for k2, v2 in list(v.items())[:3]:
            print(f"  {k2}: {type(v2).__name__} = {v2}")
    else:
        print(f"\n{k}: {type(v).__name__} = {v}")
