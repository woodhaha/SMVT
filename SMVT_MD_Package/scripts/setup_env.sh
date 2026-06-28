#!/usr/bin/env bash
# SMVT MD Environment Setup
# =========================
# One-command setup for MD simulations on Windows/Linux/Mac
# Creates conda environment with all dependencies

set -e

ENV_NAME="smvt-md"
PYTHON_VER="3.11"

echo "========================================"
echo "SMVT MD Environment Setup"
echo "========================================"

# ── 1. Check/install conda ────────────────────────────────────
if ! command -v conda &> /dev/null; then
    echo "ERROR: conda not found. Install Miniconda first:"
    echo "  https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi
echo "[✓] conda found: $(conda --version)"

# ── 2. Create environment ─────────────────────────────────────
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "[✓] Environment '${ENV_NAME}' already exists"
else
    echo "Creating conda environment '${ENV_NAME}' (Python ${PYTHON_VER})..."
    conda create -n ${ENV_NAME} python=${PYTHON_VER} -y
    echo "[✓] Environment created"
fi

# ── 3. Activate and install ───────────────────────────────────
echo ""
echo "Installing packages..."
conda run -n ${ENV_NAME} conda install -c conda-forge -y \
    openmm \
    pdbfixer \
    mdtraj \
    numpy \
    pandas \
    matplotlib \
    seaborn \
    scipy \
    jupyter \
    2>&1 | tail -5

# ── 4. pip packages ───────────────────────────────────────────
conda run -n ${ENV_NAME} pip install \
    openmmforcefields \
    openff-toolkit \
    rdkit \
    2>&1 | tail -3

# ── 5. Verify ──────────────────────────────────────────────────
echo ""
echo "Verifying installation..."
conda run -n ${ENV_NAME} python -c "
import openmm;           print(f'  openmm {openmm.__version__} OK')
import openmm.app;       print(f'  openmm.app OK')
import openmm.unit;      print(f'  openmm.unit OK')
import mdtraj;           print(f'  mdtraj {mdtraj.__version__} OK')
import numpy as np;      print(f'  numpy {np.__version__} OK')
import rdkit;            print(f'  rdkit {rdkit.__version__} OK')
try:
    from openmmforcefields.generators import SystemGenerator
    print(f'  openmmforcefields OK')
except: print(f'  openmmforcefields WARNING')
"

# ── 6. Platform-specific ──────────────────────────────────────
echo ""
case "$(uname -s)" in
    Linux*)
        echo "Platform: Linux"
        echo "GPU acceleration:"
        echo "  conda install -c conda-forge openmm-cuda"
        ;;
    Darwin*)
        echo "Platform: macOS"
        echo "GPU acceleration: Metal via OpenMM (built-in)"
        ;;
    MINGW*|MSYS*|CYGWIN*)
        echo "Platform: Windows"
        echo "  GPU not available through WSL — use native Windows conda"
        ;;
esac

echo ""
echo "========================================"
echo "Setup complete!"
echo ""
echo "Activate:  conda activate ${ENV_NAME}"
echo "Test:      python -c 'import openmm; print(openmm.Platform.getPluginLoadFailures())'"
echo "Run MD:    python scripts/run_smvt_md.py --prepare && python scripts/run_smvt_md.py --compound HYDROMORPHONE --ns 100"
echo "========================================"
