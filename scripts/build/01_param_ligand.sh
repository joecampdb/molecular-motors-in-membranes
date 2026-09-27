#!/bin/bash
# GAFF2 + AM1-BCC parameters for one motor ligand (monocation).
# Usage: 01_param_ligand.sh <NAME> <path/to/ligand.sdf>
# Output in /root/work/routeA/lig/<NAME>: lig.mol2 lig.frcmod lig.allparm.frcmod lig.lib lig.prmtop lig.inpcrd sqm.out logs
set -euo pipefail
B=/root/miniforge3/envs/md/bin
NAME=$1; SDF=$2
W=/root/work/routeA/lig/$NAME; mkdir -p "$W"; cd "$W"
cp "$SDF" input.sdf
sha256sum input.sdf > input.sdf.sha256
$B/antechamber -i input.sdf -fi mdl -o lig.mol2 -fo mol2 -c bcc -nc 1 -at gaff2 -rn LIG -s 2 > antechamber.log 2>&1
$B/parmchk2 -i lig.mol2 -f mol2 -o lig.frcmod -s gaff2 > parmchk2.log 2>&1
$B/parmchk2 -i lig.mol2 -f mol2 -o lig.allparm.frcmod -s gaff2 -a Y >> parmchk2.log 2>&1
cat > tleap.in <<T
source leaprc.gaff2
loadamberparams lig.frcmod
LIG = loadmol2 lig.mol2
check LIG
saveoff LIG lig.lib
saveamberparm LIG lig.prmtop lig.inpcrd
savepdb LIG lig_tleap.pdb
quit
T
$B/tleap -f tleap.in > tleap.log 2>&1
$B/python - <<PY
import parmed as pmd
p = pmd.load_file('lig.prmtop','lig.inpcrd')
q = sum(a.charge for a in p.atoms)
print(f"$NAME: atoms {len(p.atoms)} net charge {q:+.4f} bonds {len(p.bonds)} angles {len(p.angles)} dihedrals {len(p.dihedrals)} types {len(set(a.type for a in p.atoms))}")
PY
echo "$NAME frcmod: $(grep -c ATTN lig.frcmod || true) ATTN lines; $(grep -ci penalty lig.frcmod || true) penalty-annotated lines"
grep -iE "warning|error" tleap.log | grep -viE "^$" | head -5 || true
echo "$NAME done"
