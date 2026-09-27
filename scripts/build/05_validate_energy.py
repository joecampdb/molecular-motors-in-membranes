"""Conversion fidelity check: single-point energy decomposition of the Amber prmtop and of the
ParmEd-written GROMACS topology, both evaluated by OpenMM on identical coordinates, plus explicit
1-4 Lennard-Jones and 1-4 Coulomb sums (the terms most likely to be mangled in a conversion).

Usage: python 05_validate_energy.py <sys.prmtop> <sys.inpcrd|rst7> <sys.top> <out.json>
The recorded LJ-14 / Coulomb-14 values are the numbers a zero-step GROMACS energy evaluation of
sys.top should reproduce (GROMACS terms "LJ-14" and "Coulomb-14").
"""
import json
import sys

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import parmed as pmd
from parmed.openmm import energy_decomposition_system

prmtop, inpcrd, top, out = sys.argv[1:5]
amber = pmd.load_file(prmtop, inpcrd)
gmx = pmd.load_file(top)
gmx.coordinates = amber.coordinates
gmx.box = amber.box
names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
plat = "CUDA" if "CUDA" in names else "CPU"


def make_system(struct):
    system = struct.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometer,
                                 constraints=None, rigidWater=False, ewaldErrorTolerance=1e-5)
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
    return system


def one_four(struct, system):
    nb = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
    X = np.array(struct.coordinates) / 10.0
    lj = coul = 0.0
    n = 0
    for k in range(nb.getNumExceptions()):
        i, j, qq, sig, eps = nb.getExceptionParameters(k)
        qq = qq.value_in_unit(u.elementary_charge ** 2)
        sig = sig.value_in_unit(u.nanometer)
        eps = eps.value_in_unit(u.kilojoules_per_mole)
        if qq == 0 and eps == 0:
            continue
        r = float(np.linalg.norm(X[i] - X[j]))
        n += 1
        lj += 4 * eps * ((sig / r) ** 12 - (sig / r) ** 6)
        coul += 138.935458 * qq / r
    return n, lj, coul


result = {"platform": plat, "terms": [], "one_four": {}}
decomp = {}
for label, struct in (("amber", amber), ("gromacs_top", gmx)):
    system = make_system(struct)
    decomp[label] = {name: float(val) for name, val in energy_decomposition_system(struct, system, platform=plat, nrg=u.kilojoules_per_mole)}
    n, lj, coul = one_four(struct, system)
    result["one_four"][label] = {"pairs": n, "LJ14_kJ_mol": round(lj, 3), "Coulomb14_kJ_mol": round(coul, 3)}

worst = 0.0
for k in sorted(set(decomp["amber"]) | set(decomp["gromacs_top"])):
    a, g = decomp["amber"].get(k, 0.0), decomp["gromacs_top"].get(k, 0.0)
    rel = abs(g - a) / max(abs(a), 1.0)
    worst = max(worst, rel)
    result["terms"].append({"term": k, "amber_kJ_mol": round(a, 3), "gromacs_top_kJ_mol": round(g, 3), "diff_kJ_mol": round(g - a, 3), "rel": rel})
ta, tg = sum(decomp["amber"].values()), sum(decomp["gromacs_top"].values())
lj_a, lj_g = result["one_four"]["amber"]["LJ14_kJ_mol"], result["one_four"]["gromacs_top"]["LJ14_kJ_mol"]
result.update({"total_amber_kJ_mol": round(ta, 3), "total_gromacs_top_kJ_mol": round(tg, 3), "total_diff_kJ_mol": round(tg - ta, 3),
               "worst_relative_term_diff": worst, "LJ14_diff_kJ_mol": round(lj_g - lj_a, 3),
               "pass": bool(worst < 1e-4 and abs(lj_g - lj_a) < 1.0)})
json.dump(result, open(out, "w"), indent=2)
for r in result["terms"]:
    print(f"{r['term']:<24} amber {r['amber_kJ_mol']:>14.3f}  gmx-top {r['gromacs_top_kJ_mol']:>14.3f}  diff {r['diff_kJ_mol']:>10.3f}")
for label in ("amber", "gromacs_top"):
    o = result["one_four"][label]
    print(f"1-4 {label:<12} pairs {o['pairs']}  LJ-14 {o['LJ14_kJ_mol']:.3f}  Coulomb-14 {o['Coulomb14_kJ_mol']:.3f}")
print(f"TOTAL amber {ta:.3f} gmx-top {tg:.3f} diff {tg - ta:.3f}  worst rel {worst:.2e}  pass={result['pass']}")
