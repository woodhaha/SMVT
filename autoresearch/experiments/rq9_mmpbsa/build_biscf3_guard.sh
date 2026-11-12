#!/bin/bash
MD=/usr/local/miniconda3/envs/md_env/bin/python
GEN=/root/md/prod_mem_generic.py
L=BISCF3
WD=/root/md/build/$L
OUT=$WD/traj
CHK=$WD/prod.chk
mkdir -p "$OUT"
log() { echo "[$(date '+%F %T')] $*" >> "$WD/guard.log"; }
log "guard start"
while true; do
  timeout -k 300 54000 "$MD" "$GEN" "$L" "$WD" "$OUT" "$CHK" 100 >> "$WD/prod.log" 2>&1
  if grep -q "DONE $L" "$WD/prod.log"; then log "COMPLETED"; break; fi
  log "exited/crashed, restart in 30s"; sleep 30
done
