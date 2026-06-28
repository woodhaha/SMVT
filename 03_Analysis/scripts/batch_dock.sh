#!/bin/bash
# Batch docking: FDA drugs vs SMVT receptor
# AutoDock Vina batch runner

DATA="D:/Researching/SMVT/03_Analysis/data"
RECEPTOR="$DATA/receptor.pdbqt"
LIGANDS_DIR="$DATA/pdbqt_ligands"
RESULTS_DIR="$DATA/results"
CONFIG="$DATA/docking_config.json"

mkdir -p "$RESULTS_DIR"

# Extract config params
CX=5.03; CY=1.5; CZ=-2.44
SX=25.0; SY=25.0; SZ=25.0
EXHAUST=16
NUM_MODES=9
ERANGE=3

total=0
success=0
failed=0

echo "========================================="
echo "FDA Drug Docking: SMVT Receptor"
echo "Receptor: $(basename $RECEPTOR)"
echo "Grid: center=($CX, $CY, $CZ) size=($SX, $SY, $SZ)"
echo "Exhaustiveness: $EXHAUST | Modes: $NUM_MODES"
echo "========================================="
echo ""

for lig in "$LIGANDS_DIR"/*.pdbqt; do
    [ -f "$lig" ] || continue
    name=$(basename "$lig" .pdbqt)
    out="$RESULTS_DIR/${name}_out.pdbqt"
    log="$RESULTS_DIR/${name}.log"

    total=$((total + 1))
    echo "[$total] $name ..."

    vina --receptor "$RECEPTOR" \
         --ligand "$lig" \
         --out "$out" \
         --center_x $CX --center_y $CY --center_z $CZ \
         --size_x $SX --size_y $SY --size_z $SZ \
         --exhaustiveness $EXHAUST \
         --num_modes $NUM_MODES \
         --energy_range $ERANGE \
         > "$log" 2>&1

    if [ $? -eq 0 ]; then
        # Extract best affinity
        best=$(grep -m1 "^\s*1\s" "$log" | awk '{print $2}')
        if [ -n "$best" ]; then
            printf "    Best affinity: %.1f kcal/mol\n" "$best"
        fi
        success=$((success + 1))
    else
        echo "    FAILED"
        failed=$((failed + 1))
    fi
done

echo ""
echo "========================================="
echo "DONE: $success/$total docked ($failed failed)"
echo "Results: $RESULTS_DIR/"
echo "========================================="
