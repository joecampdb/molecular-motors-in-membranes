#!/bin/bash
# Per-replica observables with GROMACS tools only, for all production replicas of MM1 and MM2.
# For each replica: reduce the trajectory to the ligand + phosphorus atoms (molecules whole), then
#   tilt      = |90 - angle(axle C7->C14, z)|  via gmx gangle            [deg]
#   axle_z    = z(axle midpoint) - z(COM of all P)  via gmx traj -com   [A]
#   nplus_z   = z(ammonium N)   - z(COM of all P)                        [A]
#   pp        = mean z(upper-leaflet P) - mean z(lower-leaflet P)         [nm]
#   apl       = Box-X * Box-Y / (lipids per leaflet)                      [A^2]  (from the energy file)
# Output per replica: /root/work/routeA/sys/<M>/prod/r<N>/obs/series.dat  (t_ps tilt axle_z nplus_z pp)
#                     and box.xvg. Usage: 12_observables_all_gmx.sh [lipids_per_leaflet=44.5]
unset OMP_NUM_THREADS; source /root/miniforge3/envs/martini/bin/GMXRC 2>/dev/null; GMX=/root/miniforge3/envs/martini/bin/gmx
NLEAF=${1:-44.5}
one () {
  local M=$1 R=$2 D=/root/work/routeA/sys/$1/prod/r$2
  [ -f $D/prod.gro ] || { echo "$M r$R: not finished, skipped"; return; }
  mkdir -p $D/obs && cd $D/obs || return
  $GMX select -s ../prod.tpr -on small.ndx -select '"SMALL" resname LIG or name P31' > /dev/null 2>&1
  echo SMALL | $GMX trjconv -f ../prod.xtc -s ../prod.tpr -n small.ndx -pbc mol -o small.xtc > trjconv.log 2>&1
  echo SMALL | $GMX convert-tpr -s ../prod.tpr -n small.ndx -o small.tpr > /dev/null 2>&1
  $GMX select -s small.tpr -on obs.ndx -select '"AXLE" atomnr 7 15' '"PHOS" name P31' '"NPLUS" atomnr 32' > /dev/null 2>&1
  $GMX gangle -f small.xtc -s small.tpr -n obs.ndx -g1 vector -group1 AXLE -g2 z -oav tilt.xvg > gangle.log 2>&1
  echo AXLE  | $GMX traj -f small.xtc -s small.tpr -n obs.ndx -com -ox axle_com.xvg  > /dev/null 2>&1
  echo NPLUS | $GMX traj -f small.xtc -s small.tpr -n obs.ndx -com -ox nplus_com.xvg > /dev/null 2>&1
  echo PHOS  | $GMX traj -f small.xtc -s small.tpr -n obs.ndx -com -ox phos_com.xvg  > /dev/null 2>&1
  echo PHOS  | $GMX traj -f small.xtc -s small.tpr -n obs.ndx -ox phos_all.xvg       > /dev/null 2>&1
  printf "Box-X\nBox-Y\n" | $GMX energy -f ../prod.edr -o box.xvg > /dev/null 2>&1
  paste <(grep -vE "^[#@]" tilt.xvg | awk '{a=$2; t=a>90?a-90:90-a; print $1, t}') \
        <(grep -vE "^[#@]" axle_com.xvg | awk '{print $4}') \
        <(grep -vE "^[#@]" nplus_com.xvg | awk '{print $4}') \
        <(grep -vE "^[#@]" phos_com.xvg | awk '{print $4}') \
        <(grep -vE "^[#@]" phos_all.xvg | awk '{n=0; for(i=4;i<=NF;i+=3){z[n++]=$i}; s=0; for(i=0;i<n;i++)s+=z[i]; m=s/n; su=0;nu=0;sl=0;nl=0; for(i=0;i<n;i++){if(z[i]>m){su+=z[i];nu++}else{sl+=z[i];nl++}}; print su/nu-sl/nl}') \
    | awk '{printf "%.0f %.3f %.3f %.3f %.4f\n", $1, $2, ($3-$5)*10, ($4-$5)*10, $6}' > series.dat
  rm -f small.xtc
  echo "$M r$R: $(wc -l < series.dat) frames, last $(tail -1 series.dat | awk '{print $1/1000}') ns"
}
for M in MM1 MM2; do for R in 1 2 3 4 5; do one $M $R & done; done; wait
echo "all replicas processed"
