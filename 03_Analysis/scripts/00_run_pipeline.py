#!/usr/bin/env python3
"""
00_run_pipeline.py
Master pipeline for SMVT (SLC5A6) virtual screening workflow.

Run this script to execute the complete pipeline:
1. Prepare known interactor and FDA-approved drug ligands
2. Analyze SMVT structure and define docking pocket
3. Convert receptor PDB to PDBQT
4. Run AutoDock Vina virtual screening
5. Analyze results and generate report
"""

import os, sys, subprocess, time
from pathlib import Path

PYTHON = r"C:\anaconda3\python.exe"
SCRIPTS_DIR = Path(r"D:\Researching\SMVT\03_Analysis\scripts")

def run_step(name, script):
    """Run a Python script step."""
    print(f"\n{'='*60}")
    print(f"STEP: {name}")
    print(f"{'='*60}\n")

    script_path = SCRIPTS_DIR / script
    if not script_path.exists():
        print(f"ERROR: Script not found: {script_path}")
        return False

    start = time.time()
    try:
        result = subprocess.run(
            [PYTHON, str(script_path)],
            capture_output=True, text=True, timeout=3600
        )
        elapsed = time.time() - start

        print(result.stdout)
        if result.stderr:
            print(f"STDERR: {result.stderr[:500]}")

        if result.returncode == 0:
            print(f"\n  Completed in {elapsed:.0f}s")
            return True
        else:
            print(f"\n  FAILED (exit={result.returncode}) in {elapsed:.0f}s")
            return False
    except subprocess.TimeoutExpired:
        print(f"\n  TIMEOUT after {(time.time()-start):.0f}s")
        return False
    except Exception as e:
        print(f"\n  ERROR: {e}")
        return False


def main():
    print("=" * 60)
    print("SMVT (SLC5A6) Virtual Screening Pipeline")
    print("=" * 60)

    pipeline = [
        ("Prepare ligands", "01_prepare_ligands.py"),
        ("Analyze structure", "02_analyze_structure.py"),
        ("Convert receptor to PDBQT", "pdb_to_pdbqt.py"),
    ]

    for name, script in pipeline:
        if not run_step(name, script):
            print(f"\nPipeline ABORTED at step: {name}")
            sys.exit(1)

    # Step 4: Run docking (can be long)
    print(f"\n{'='*60}")
    print(f"FINAL STEP: Virtual Screening Docking")
    print(f"{'='*60}")
    print(f"Running 03_run_docking.py...")

    start = time.time()
    result = subprocess.run(
        [PYTHON, str(SCRIPTS_DIR / "03_run_docking.py")],
        capture_output=True, text=True, timeout=86400  # 24h timeout for all docking
    )
    elapsed = time.time() - start

    print(result.stdout[-3000:] if len(result.stdout) > 3000 else result.stdout)
    if result.stderr:
        print(f"STDERR (last 1000): {result.stderr[-1000:]}")

    if result.returncode == 0:
        print(f"\n{'='*60}")
        print(f"PIPELINE COMPLETE! ({elapsed:.0f}s total)")
        print(f"{'='*60}")

        # Show output files
        outputs = Path(r"D:\Researching\SMVT\03_Analysis\outputs")
        print("\nOutput files:")
        for f in outputs.glob("*"):
            print(f"  {f.name} ({os.path.getsize(f):,} bytes)")
    else:
        print(f"\nPIPELINE FAILED at docking step (exit={result.returncode})")
        print(f"Time elapsed: {elapsed:.0f}s")


if __name__ == "__main__":
    main()
