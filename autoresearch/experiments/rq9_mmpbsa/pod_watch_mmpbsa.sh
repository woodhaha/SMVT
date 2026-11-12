#!/bin/bash
# Pod-side watcher: poll for BISCF3 MD completion, then auto-run MM-GBSA.
# Runs detached on the pod; independent of any local machine.
WD=/root/md/build/BISCF3
LOG=$WD/rq9_watch.log
echo "[$(date '+%F %T')] watcher start" >> "$LOG"
while true; do
  if grep -q "DONE BISCF3" "$WD/prod.log" 2>/dev/null; then
    echo "[$(date '+%F %T')] MD done -> running MM-GBSA" >> "$LOG"
    bash /root/md/rq9/run_mmpbsa.sh >> "$LOG" 2>&1
    echo "[$(date '+%F %T')] MM-GBSA finished" >> "$LOG"
    break
  fi
  sleep 900
done
