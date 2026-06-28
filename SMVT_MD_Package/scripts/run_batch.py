#!/usr/bin/env python3
"""
SMVT MD Batch Runner — Top 3 + 3 Controls, 50ns each, GPU sequential
====================================================================
Single GPU → one simulation at a time. 6 × 50ns = 300ns total.
RTX 2000 Ada: ~40ns/day (OpenCL) → ~7.5 days | ~70ns/day (CUDA) → ~4 days
"""

import os, sys, json, time, subprocess, logging
from datetime import datetime, timedelta

PACKAGE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOG_FILE = os.path.join(PACKAGE_DIR, "batch_md.log")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.FileHandler(LOG_FILE), logging.StreamHandler(sys.stdout)]
)
log = logging.getLogger(__name__)

# ── Run order: top 3 + 3 controls ─────────────────────────────
COMPOUNDS = [
    ("HYDROMORPHONE",       50, "Top 1 — Opioid"),
    ("FUROSEMIDE",           50, "Top 2 — Diuretic"),
    ("NAFTAZONE",            50, "Top 3 — Hemostatic"),
    ("BIOTIN",               50, "Reference — Natural substrate"),
    ("GABAPENTIN_ENACARBIL", 50, "Positive control — FDA SMVT drug"),
    ("RIBOFLAVIN",           50, "Negative control — Non-binding vitamin"),
]

# Use working conda env (OpenCL GPU, openmm 8.5.2 with all force fields)
PYTHON = r"C:\anaconda3\envs\smvt-md\python.exe"
if not os.path.exists(PYTHON):
    PYTHON = os.path.join(PACKAGE_DIR, ".venv", "Scripts", "python.exe")
if not os.path.exists(PYTHON):
    PYTHON = sys.executable  # last resort
SCRIPT = os.path.join(PACKAGE_DIR, "scripts", "run_smvt_md.py")


def detect_platform():
    """Check which OpenMM platform is available."""
    import openmm as mm
    for name in ["CUDA", "OpenCL", "CPU"]:
        try:
            mm.Platform.getPlatformByName(name)
            return name
        except:
            pass
    return "Reference"


def run_one(compound_id, ns):
    """Run MD for one compound, return elapsed time."""
    log.info(f"{'='*60}")
    log.info(f"START: {compound_id} ({ns}ns) on {detect_platform()}")
    log.info(f"{'='*60}")

    t0 = time.time()
    result = subprocess.run(
        [PYTHON, SCRIPT, "--compound", compound_id, "--ns", str(ns)],
        cwd=PACKAGE_DIR,
        capture_output=False,
        text=True,
    )
    elapsed = time.time() - t0

    if result.returncode == 0:
        log.info(f"DONE: {compound_id} in {elapsed/3600:.1f}h ({ns/(elapsed/86400):.1f} ns/day)")
        return True, elapsed
    else:
        log.error(f"FAILED: {compound_id} (exit code {result.returncode})")
        return False, elapsed


def main():
    platform = detect_platform()
    log.info(f"Platform: {platform}")
    log.info(f"Compounds: {len(COMPOUNDS)} | Total: {sum(n for _,n,_ in COMPOUNDS)}ns")
    log.info(f"Estimated: {sum(n for _,n,_ in COMPOUNDS)/40 if platform=='OpenCL' else sum(n for _,n,_ in COMPOUNDS)/70:.0f} days ({platform})")

    results = []
    total_start = time.time()

    for i, (cid, ns, desc) in enumerate(COMPOUNDS):
        log.info(f"\n[{i+1}/6] {cid} — {desc}")
        ok, elapsed = run_one(cid, ns)
        results.append({
            "compound": cid, "ns": ns, "success": ok,
            "elapsed_h": elapsed/3600, "ns_per_day": ns/(elapsed/86400) if elapsed > 0 else 0,
            "platform": platform,
        })

        # Save progress
        with open(os.path.join(PACKAGE_DIR, "batch_progress.json"), "w") as f:
            json.dump({"completed": i+1, "total": 6, "results": results,
                       "started": datetime.now().isoformat()}, f, indent=2)

    total_elapsed = (time.time() - total_start) / 3600
    total_ns = sum(r["ns"] for r in results)
    log.info(f"\n{'='*60}")
    log.info(f"BATCH COMPLETE: {total_elapsed:.1f}h for {total_ns}ns")
    log.info(f"Average speed: {total_ns/(total_elapsed*24):.1f} ns/day")
    log.info(f"{'='*60}")

    # Print summary
    print("\n" + "="*60)
    print("SMVT MD BATCH RESULTS")
    print("="*60)
    for r in results:
        status = "✅" if r["success"] else "❌"
        print(f"  {status} {r['compound']:<25} {r['ns']}ns  {r['elapsed_h']:.1f}h  {r['ns_per_day']:.0f}ns/day")
    print(f"\nTotal: {total_elapsed:.1f}h ({total_elapsed/24:.1f} days)")


if __name__ == "__main__":
    main()
