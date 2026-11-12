#!/bin/bash
# Run MM-GBSA on COOH-bisCF3 vs SMVT (single trajectory, GB-OBC2).
# Run AFTER production MD completes (traj/BISCF3_100ns.dcd exists).
set -e
MD=/usr/local/miniconda3/envs/md_env/bin
WD=/root/md/build/BISCF3
cd "$WD"

# 1. dry topologies (strip water/lipids/ions; keep protein + ligand)
$MD/ante-MMPBSA.py -p system_full.prmtop \
  -c complex_dry.prmtop -r receptor_dry.prmtop -l ligand_dry.prmtop \
  -s ":WAT,Na+,Cl-,PA,OL,PC" -m ":UNL" --radii=mbondi2

# 2. MM-GBSA
$MD/MMPBSA.py -O -i /root/md/rq9/mmpbsa.in -o FINAL_RESULTS_MMPBSA.dat \
  -sp system_full.prmtop -cp complex_dry.prmtop -rp receptor_dry.prmtop -lp ligand_dry.prmtop \
  -y traj/BISCF3_100ns.dcd

echo "=== MM-GBSA done ==="
grep -A6 "DELTA TOTAL" FINAL_RESULTS_MMPBSA.dat | head -20
