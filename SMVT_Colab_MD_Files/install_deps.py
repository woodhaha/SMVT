"""Install MD dependencies on Colab runtime."""
import subprocess, sys, importlib

PACKAGES = [
    "openmm",
    "pdbfixer",
    "rdkit",
    "mdtraj",
    "numpy",
    "matplotlib",
    "openmmforcefields",
]

print("=" * 50)
print("Installing MD dependencies on Colab...")
print("=" * 50)

for pkg in PACKAGES:
    try:
        importlib.import_module(pkg.replace("-", "_"))
        print(f"  {pkg}: already installed")
    except ImportError:
        print(f"  {pkg}: installing...")
        subprocess.check_call(
            [sys.executable, "-m", "pip", "install", "-q", pkg],
            timeout=300
        )
        print(f"  {pkg}: done")

# Verify OpenMM GPU
import openmm as mm
print(f"\nOpenMM {mm.__version__}")
for pname in ["CUDA", "OpenCL", "CPU"]:
    try:
        plat = mm.Platform.getPlatformByName(pname)
        print(f"  GPU Platform: {plat.getName()} ✓")
        break
    except:
        continue

print("\nAll dependencies ready!")
