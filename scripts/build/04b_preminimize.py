"""Remove packing overlaps before conversion: L-BFGS minimisation of the packed Amber system with OpenMM,
in the Amber representation, so GROMACS never starts from near-coincident atoms.

Usage: python 04b_preminimize.py <prmtop> <crd> <out.rst7> <report.json>
Reports energies before/after, the closest non-bonded pair before/after, and how far atoms moved.
"""
import json
import sys

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import parmed as pmd
from scipy.spatial import cKDTree

prmtop, crd, out, rep = sys.argv[1:5]
s = pmd.load_file(prmtop, crd)
bonded = set((min(b.atom1.idx, b.atom2.idx), max(b.atom1.idx, b.atom2.idx)) for b in s.bonds)


def closest(X, cutoff=1.5):
    pairs = cKDTree(X).query_pairs(cutoff, output_type="ndarray")
    keep = [k for k, (i, j) in enumerate(pairs) if (min(i, j), max(i, j)) not in bonded]
    if not keep:
        return {"n_nonbonded_pairs_below_cutoff": 0, "min_distance_A": None}
    p = pairs[keep]
    d = np.linalg.norm(X[p[:, 0]] - X[p[:, 1]], axis=1)
    k = int(np.argmin(d))
    return {"n_nonbonded_pairs_below_cutoff": int(len(d)), "n_below_1A": int((d < 1.0).sum()),
            "min_distance_A": round(float(d[k]), 3),
            "closest_pair": [f"{s.atoms[p[k,0]].residue.name}{s.atoms[p[k,0]].residue.idx}:{s.atoms[p[k,0]].name}",
                             f"{s.atoms[p[k,1]].residue.name}{s.atoms[p[k,1]].residue.idx}:{s.atoms[p[k,1]].name}"]}


X0 = np.array(s.coordinates, dtype=float)
before = closest(X0)
system = s.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometer, constraints=None,
                        rigidWater=False, ewaldErrorTolerance=1e-5)
names = [mm.Platform.getPlatform(i).getName() for i in range(mm.Platform.getNumPlatforms())]
plat = mm.Platform.getPlatformByName("CUDA" if "CUDA" in names else "CPU")
props = {"Precision": "mixed"} if plat.getName() == "CUDA" else {}
ctx = mm.Context(system, mm.VerletIntegrator(0.001 * u.picoseconds), plat, props)
ctx.setPeriodicBoxVectors(*s.box_vectors)
ctx.setPositions(s.positions)
e0 = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
mm.LocalEnergyMinimizer.minimize(ctx, tolerance=10.0 * u.kilojoules_per_mole / u.nanometer, maxIterations=10000)
st = ctx.getState(getEnergy=True, getPositions=True)
e1 = st.getPotentialEnergy().value_in_unit(u.kilojoules_per_mole)
X1 = st.getPositions(asNumpy=True).value_in_unit(u.angstrom)
after = closest(X1)
disp = np.linalg.norm(X1 - X0, axis=1)
heavy = np.array([a.atomic_number > 1 for a in s.atoms])
s.coordinates = X1
s.save(out, overwrite=True)
report = {"prmtop": prmtop, "input_coordinates": crd, "output_coordinates": out, "platform": plat.getName(),
          "potential_kJ_mol_before": e0, "potential_kJ_mol_after": e1,
          "closest_nonbonded_before": before, "closest_nonbonded_after": after,
          "heavy_atom_displacement_A": {"rmsd": round(float(np.sqrt((disp[heavy] ** 2).mean())), 3),
                                        "max": round(float(disp[heavy].max()), 3)},
          "note": "Unrestrained L-BFGS minimisation in the Amber representation; removes packmol overlaps only."}
json.dump(report, open(rep, "w"), indent=2)
print(json.dumps(report))
