"""Write the xtb geometry and scan inputs for the MM1 cation.

Atom order comes from lig.prmtop, so xtb and the GAFF2 force field index atoms identically and the
two scans can be compared directly.

Two torsions:
  twist   C3-C7=C14-C15, the overcrowded alkene itself. Both alkene carbons are ring-fused, so this
          coordinate twists the whole molecule and sets its helicity. It carries the parameters
          parmchk2 filled by analogy with the highest penalties, so it is the one worth testing.
  arylN   the rotatable aryl-N bond of the arylamine, which orients the tail that carries the
          ammonium anchor into the head-group region.
"""
import os
import sys

import numpy as np
import parmed as pmd

out = sys.argv[1]
os.makedirs(out, exist_ok=True)
p = pmd.load_file("/root/work/routeA/lig/MM1/lig.prmtop", "/root/work/routeA/lig/MM1/lig.inpcrd")
names = [a.name for a in p.atoms]
elem = [a.element_name for a in p.atoms]
xyz = np.array(p.coordinates)
nb = {i: [] for i in range(len(p.atoms))}
for b in p.bonds:
    nb[b.atom1.idx].append(b.atom2.idx)
    nb[b.atom2.idx].append(b.atom1.idx)
heavy = lambda i: [j for j in nb[i] if elem[j] != "H"]
ring = {a.idx for a in p.atoms if len(heavy(a.idx)) >= 2 and elem[a.idx] == "C"}


def dihedral(q, c=xyz):
    a, b, cc, d = [c[i] for i in q]
    b0, b1, b2 = a - b, cc - b, d - cc
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


c7, c14 = names.index("C7"), names.index("C14")
twist = (heavy(c7)[0] if heavy(c7)[0] != c14 else heavy(c7)[1], c7, c14,
         heavy(c14)[0] if heavy(c14)[0] != c7 else heavy(c14)[1])

# arylamine nitrogen: the neutral N, bonded to one aromatic carbon and one sp3 carbon
nitrogens = [i for i, e in enumerate(elem) if e == "N"]
aryl_n = None
for n in nitrogens:
    hv = heavy(n)
    if len(hv) == 2:                       # N1: aryl + CH2  (the ammonium N has three carbons)
        car = [j for j in hv if len(heavy(j)) >= 3 and j in ring]
        if car:
            aryl_n = (car[0], n, [j for j in hv if j != car[0]][0])
            break
assert aryl_n, "arylamine nitrogen not found"
car, n1, ch2 = aryl_n
ring_c = [j for j in heavy(car) if j != n1][0]
arylN = (ring_c, car, n1, ch2)

scans = {"twist": {"atoms": twist, "start": -70.0, "end": 70.0, "steps": 29},
         "arylN": {"atoms": arylN, "start": -180.0, "end": 180.0, "steps": 37}}

with open(os.path.join(out, "mm1.xyz"), "w") as f:
    f.write("%d\nMM1 cation, atom order from lig.prmtop\n" % len(p.atoms))
    for e, c in zip(elem, xyz):
        f.write("%-2s %14.8f %14.8f %14.8f\n" % (e, c[0], c[1], c[2]))

for key, s in scans.items():
    q = s["atoms"]
    s["names"] = [names[i] for i in q]
    s["current_deg"] = dihedral(q)
    with open(os.path.join(out, "scan_%s.inp" % key), "w") as f:
        f.write("$constrain\n  force constant=1.5\n  dihedral: %d,%d,%d,%d,auto\n" % tuple(i + 1 for i in q))
        f.write("$scan\n  mode=concerted\n  1: %.1f,%.1f,%d\n$end\n" % (s["start"], s["end"], s["steps"]))
    print("%-6s %-22s now %+7.2f deg   scan %+.0f..%+.0f in %d steps"
          % (key, "-".join(s["names"]), s["current_deg"], s["start"], s["end"], s["steps"]))

np.save(os.path.join(out, "scans.npy"), np.array(
    [[key] + [str(i) for i in s["atoms"]] for key, s in scans.items()], dtype=object), allow_pickle=True)
print("net charge %.3f -> xtb --chrg 1" % sum(a.charge for a in p.atoms))
