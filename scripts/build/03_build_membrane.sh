#!/bin/bash
# Pack one motor ligand into a POPE:POPG (1:3) bilayer with TIP3P water, Na+ neutralisation and
# 150 mM NaCl, then parametrise with tleap (Lipid21 + GAFF2 ligand + TIP3P) via packmol-memgen.
# Usage: 03_build_membrane.sh <NAME>
# Expects /root/work/routeA/lig/<NAME>/{lig.frcmod,lig.lib} and .../lig_oriented.pdb
# Output: /root/work/routeA/sys/<NAME>/ with packed PDB, prmtop/inpcrd and logs.
set -euo pipefail
B=/root/miniforge3/envs/md/bin
export AMBERHOME=/root/miniforge3/envs/md
export PATH=$B:$PATH
NAME=$1
L=/root/work/routeA/lig/$NAME
W=/root/work/routeA/sys/$NAME; mkdir -p "$W"; cd "$W"
cp "$L/lig.frcmod" "$L/lig.lib" "$L/lig_oriented.pdb" .
# Lateral patch 53 A gives about 44 lipids per leaflet at the areas per lipid used by packmol-memgen
# (target from the source paper: 11 POPE + 33 POPG per leaflet). Water layer 50 A each side of the bilayer.
# packmol-memgen's --saltcon is the TOTAL cation concentration. Neutralising the POPG charge alone already
# exceeds 0.15 M, so pass 1 only reads the neutralisation requirement and pass 2 asks for that plus 0.15 M
# NaCl, which reproduces "Na+ neutralisation and 150 mM NaCl" in the tool's own volume convention.
COMMON=(--pdb lig_oriented.pdb --preoriented --notprotonate --nottrim --charge_pdb_delta 1
        --lipids POPE:POPG --ratio 1:3 --distxy_fix 53 --dist_wat 50
        --salt --salt_c Na+ --salt_a Cl-
        --parametrize --ffwat tip3p --fflip lipid21 --ffprot ff14SB
        --ligand_param lig.frcmod:lig.lib --gaff2
        --nloop 40 --noprogress --overwrite -o packed_${NAME}.pdb)
$B/packmol-memgen "${COMMON[@]}" --saltcon 0.15 --log pass1.log > pass1.stdout 2>&1 || true
NEUT=$(grep -m1 "Positive ion concentration" pass1.stdout | awk '{print $NF}')
if [ -z "$NEUT" ]; then echo "pass 1 did not report the neutralisation concentration"; tail -20 pass1.stdout; exit 1; fi
SALTCON=$(python3 -c "print(round($NEUT + 0.15, 3))")
echo "neutralisation requires $NEUT M cations; requesting total $SALTCON M (neutralisation + 0.15 M NaCl)"
$B/packmol-memgen "${COMMON[@]}" --saltcon "$SALTCON" --log packmol-memgen.log > build.stdout 2>&1
ls -la *.prmtop *.inpcrd *.pdb 2>/dev/null
echo "lipid/water/ion summary from tleap output:"
grep -E "Total (atoms|residues)|Added|Warning|WARNING|Error" -i leap.log 2>/dev/null | head -12 || true
echo "$NAME build done"
