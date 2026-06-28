#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# AutoDL RTX 4090 Setup — SMVT MD Simulation
# Run once after connecting to AutoDL instance
# ═══════════════════════════════════════════════════════════════
set -e

echo "=== SMVT MD — AutoDL RTX 4090 Setup ==="

# 1. Academic network acceleration
source /etc/network_turbo 2>/dev/null && echo "[OK] Network turbo enabled" || echo "[SKIP] No network_turbo"

# 2. Create conda environment (if not exists)
if conda env list | grep -q smvt-md; then
    echo "[OK] conda env smvt-md already exists"
else
    echo "[...] Creating conda environment (this takes ~5min)..."
    conda create -n smvt-md python=3.11 -y
    source activate smvt-md
    conda install -c conda-forge openmm openff-toolkit pdbfixer rdkit numpy matplotlib seaborn -y
    pip install openmmforcefields mdtraj
    echo "[OK] Environment created"
fi

# 3. Verify GPU
echo ""
echo "=== GPU Info ==="
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

# 4. Verify OpenMM sees the GPU
source activate smvt-md
python -c "
import openmm as mm
for i in range(mm.Platform.getNumPlatforms()):
    p = mm.Platform.getPlatform(i)
    print(f'Platform {i}: {p.getName()}')
print('GPU OK' if any('CUDA' in mm.Platform.getPlatform(i).getName() or 'OpenCL' in mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())) else 'GPU NOT FOUND!')
"

echo ""
echo "=== Setup Complete ==="
echo "Next: bash run_hydromorphone_100ns.sh"
