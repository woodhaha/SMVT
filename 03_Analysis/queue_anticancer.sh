#!/bin/bash
# Queue: 10569 run -> anticancer docking -> merge. Detached via setsid -f.
set -u
CK=/root/docking/results.csv.checkpoint.json
cd /root/docking || exit 1

echo "$(date): waiting for combined run >= 10569"
while [ "$(python3 -c "import json;c=json.load(open('$CK'));print(len(c.get('completed',[]))+len(c.get('failed',[])))" 2>/dev/null || echo 0)" -lt 10569 ]; do
    sleep 300
done

echo "$(date): combined done. docking anticancer library (2474)"
python3 dock.py --input anticancer_new_input.csv --output anticancer_results.csv >> anticancer.log 2>&1

echo "$(date): anticancer docking done. merging"
python3 merge_anticancer.py >> anticancer.log 2>&1
echo "$(date): ALL DONE"
