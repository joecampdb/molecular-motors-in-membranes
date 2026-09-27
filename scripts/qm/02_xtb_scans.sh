#!/bin/bash
cd /root/work/qm || exit 1
export OMP_NUM_THREADS=8
export XTBPATH=/root/mamba/envs/xtb/share/xtb
XTB=/root/mamba/envs/xtb/bin/xtb
for s in twist arylN; do
  echo "[$(date +%H:%M:%S)] xtb relaxed scan: $s"
  rm -rf "$s"; mkdir -p "$s"; cd "$s" || exit 1
  $XTB ../mm1.xyz --opt --input "../scan_${s}.inp" --chrg 1 --gfn 2 --parallel 8 > xtb.log 2>&1
  rc=$?
  n=$(grep -c "^ *65 *$" xtbscan.log 2>/dev/null || echo 0)
  echo "  exit=$rc  frames=$n"
  grep -iE "convergence|failed|error" xtb.log | tail -2
  cd ..
done
echo "[$(date +%H:%M:%S)] scans done"
