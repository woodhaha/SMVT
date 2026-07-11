#!/bin/bash
# Run APBS for all 8 compounds in WSL
set -e
export PATH=/tmp/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
MM="micromamba run -n mmgbsa"
DIR=/home/woodhaha00/smvt_data/pqr
RESULTS=/home/woodhaha00/smvt_data/pqr

COMPS=("BIOTIN" "HYDROMORPHONE" "GABAPENTIN_ENACARBIL" "NAFTAZONE"
       "ESKETAMINE" "FUROSEMIDE" "PHENOBARBITAL" "RIBOFLAVIN")

for COMPOUND in "${COMPS[@]}"; do
  echo "=== $COMPOUND ==="

  # Run APBS for complex, receptor, ligand
  for SUFFIX in complex receptor ligand; do
    PQR=${DIR}/${COMPOUND}_${SUFFIX}.pqr
    IN=${DIR}/${COMPOUND}_${SUFFIX}.apbs
    OUT=${DIR}/${COMPOUND}_${SUFFIX}_apbs.out

    cat > $IN << EOF
read
  mol pqr $PQR
end
elec name ${SUFFIX}
  mg-manual
  dime 129 129 129
  glen 80 80 80
  gcent mol 1
  mol 1
  lpbe
  bcfl mdh
  sdens 10.0
  pdie 1.0
  sdie 80.0
  srad 1.4
  swin 0.3
  temp 300.0
  ion charge 1.0 conc 0.150 radius 2.0
  ion charge -1.0 conc 0.150 radius 2.0
  calcenergy total
  calcforce no
end
quit
EOF

    $MM apbs $IN > $OUT 2>&1
    echo "  $SUFFIX done"
  done

  # Parse energies
  G_COMPLEX=$(grep "Global net" ${DIR}/${COMPOUND}_complex_apbs.out 2>/dev/null | awk '{print $NF}')
  G_RECEPTOR=$(grep "Global net" ${DIR}/${COMPOUND}_receptor_apbs.out 2>/dev/null | awk '{print $NF}')
  G_LIGAND=$(grep "Global net" ${DIR}/${COMPOUND}_ligand_apbs.out 2>/dev/null | awk '{print $NF}')

  echo "  G_complex=$G_COMPLEX  G_receptor=$G_RECEPTOR  G_ligand=$G_LIGAND"

  # Handle missing values
  G_COMPLEX=${G_COMPLEX:-0}
  G_RECEPTOR=${G_RECEPTOR:-0}
  G_LIGAND=${G_LIGAND:-0}

  DG_KB=$(echo "scale=4; $G_COMPLEX - $G_RECEPTOR - $G_LIGAND" | bc -l)
  DG_KCAL=$(echo "scale=4; $DG_KB / 4.184" | bc -l)

  echo "$COMPOUND: ΔG_PBSA = ${DG_KB} kJ = ${DG_KCAL} kcal/mol"

  # Save JSON
  cat > ${DIR}/${COMPOUND}_mmpbsa.json << EOF
{
  "compound": "$COMPOUND",
  "G_complex_kJ": $G_COMPLEX,
  "G_receptor_kJ": $G_RECEPTOR,
  "G_ligand_kJ": $G_LIGAND,
  "dG_kJ": $(printf "%.2f" $DG_KB),
  "dG_kcal": $(printf "%.2f" $DG_KCAL)
}
EOF
  echo "  -> ${COMPOUND}_mmpbsa.json"
done

# Final summary
echo ""
echo "=========================================="
echo "  MM-PBSA Summary"
echo "=========================================="
printf "  %-25s %12s %12s\n" "Compound" "ΔG(kJ)" "ΔG(kcal)"
printf "  %-25s %12s %12s\n" "-------------------------" "----------" "----------"
for COMPOUND in "${COMPS[@]}"; do
  if [ -f ${DIR}/${COMPOUND}_mmpbsa.json ]; then
    DG_KJ=$(grep "dG_kJ" ${DIR}/${COMPOUND}_mmpbsa.json | grep -oP '\d+\.?\d*')
    DG_KCAL=$(grep "dG_kcal" ${DIR}/${COMPOUND}_mmpbsa.json | grep -oP '\d+\.?\d*')
    printf "  %-25s %10.2f kJ %10.2f kcal\n" "$COMPOUND" "$DG_KJ" "$DG_KCAL"
  fi
done
echo "=========================================="
