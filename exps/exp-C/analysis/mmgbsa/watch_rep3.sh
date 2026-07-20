#!/bin/bash
# Watch for exp-C rep3 completion, then extract dry trajectories and launch rep3 MM-GBSA runs.
set -u
ROOT=/home/scroll/zzhang/semaglutide-in-silico
BASE=$ROOT/exps/exp-C/analysis/mmgbsa
export PATH=/home/scroll/miniforge3/envs/cgas-md/bin:$PATH

echo "[watcher] waiting for both rep3 md_output.log to say Done..."
while true; do
  m=$(grep -c "^Done" $ROOT/exps/exp-C/md/c18_monoacid/rep3/md_output.log 2>/dev/null || true)
  d=$(grep -c "^Done" $ROOT/exps/exp-C/md/c18_diacid/rep3/md_output.log 2>/dev/null || true)
  if [ "${m:-0}" -ge 1 ] && [ "${d:-0}" -ge 1 ]; then break; fi
  sleep 300
done
echo "[watcher] $(date) rep3 MD done; extracting dry trajectories"
for s in c18_monoacid c18_diacid; do
  (cd $BASE/$s/rep3 && cpptraj -i extract.cpptraj > extract.log 2>&1)
  echo "[watcher] $s rep3 extract: $(grep -oE 'reading [0-9]+ of [0-9]+' $BASE/$s/rep3/extract.log | head -1)"
done
echo "[watcher] $(date) launching rep3 MM-GBSA runs"
cd $BASE
for s in c18_monoacid c18_diacid; do
  setsid nohup ./run_rep.sh $s 3 14 > /dev/null 2>&1 < /dev/null &
done
echo "[watcher] rep3 runs launched"
