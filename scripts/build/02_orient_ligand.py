"""Place one motor ligand in an explicit, documented starting pose for packmol-memgen.

Pose (an imposed initial condition, not a prediction):
  * axle C=C (the only non-aromatic C=C) along +x, i.e. parallel to the membrane plane;
  * protonated dimethylammonium N on the +z side of the axle (toward the upper leaflet head groups);
  * that N placed at (0, 0, +Z_N) with Z_N = 17 A, near the upper phosphate/glycerol region of a
    POPE/POPG bilayer whose midplane is z = 0.
Coordinates come from the source SDF (same atom order as the antechamber mol2); atom and residue
names come from the mol2 so tleap can match the ligand library.

Usage: python 02_orient_ligand.py <NAME> <ligand.sdf> <lig.mol2> <out.pdb> [Z_N]
"""
import json
import sys

import numpy as np
import parmed as pmd
from rdkit import Chem

name, sdf, mol2, out_pdb = sys.argv[1:5]
z_n = float(sys.argv[5]) if len(sys.argv) > 5 else 17.0

mol = Chem.MolFromMolFile(sdf, removeHs=False)
X = np.array(mol.GetConformer().GetPositions(), dtype=float)
axle = [b for b in mol.GetBonds()
        if b.GetBondType() == Chem.BondType.DOUBLE and not b.GetIsAromatic()
        and b.GetBeginAtom().GetSymbol() == "C" and b.GetEndAtom().GetSymbol() == "C"]
assert len(axle) == 1, f"expected one non-aromatic C=C, found {len(axle)}"
i, j = axle[0].GetBeginAtomIdx(), axle[0].GetEndAtomIdx()
nplus = [a.GetIdx() for a in mol.GetAtoms() if a.GetFormalCharge() == 1 and a.GetSymbol() == "N"]
assert len(nplus) == 1, f"expected one N+, found {nplus}"
k = nplus[0]


def rot_a_to_b(a, b):
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


mid = 0.5 * (X[i] + X[j])
Y = X - mid
R1 = rot_a_to_b(Y[j] - Y[i], np.array([1.0, 0, 0]))
Y = Y @ R1.T
v = Y[k].copy(); v[0] = 0.0                     # projection of axle->N+ onto the yz plane
R2 = rot_a_to_b(v, np.array([0, 0, 1.0]))       # rotation about x (v has no x component)
Y = Y @ R2.T
Y[:, 2] += z_n - Y[k, 2]                        # put N+ at z = z_n; axle midpoint stays at x = y = 0

axle_vec = Y[j] - Y[i]
tilt = np.degrees(np.arcsin(abs(axle_vec[2]) / np.linalg.norm(axle_vec)))
struct = pmd.load_file(mol2, structure=True)
assert len(struct.atoms) == len(Y), "mol2 and SDF atom counts differ"
for at, sym in zip(struct.atoms, [a.GetSymbol() for a in mol.GetAtoms()]):
    assert at.element_name.upper() == sym.upper(), f"atom order mismatch at {at.idx}: {at.element_name} vs {sym}"
struct.coordinates = Y
struct.save(out_pdb, overwrite=True)
report = {
    "ligand": name, "sdf": sdf, "axle_atoms_zero_based": [i, j], "ammonium_N_zero_based": k,
    "axle_tilt_from_membrane_plane_deg": round(tilt, 3), "axle_midpoint_A": [0.0, 0.0, round(float(0.5 * (Y[i, 2] + Y[j, 2])), 3)],
    "ammonium_N_position_A": [round(float(c), 3) for c in Y[k]],
    "z_extent_A": [round(float(Y[:, 2].min()), 2), round(float(Y[:, 2].max()), 2)],
    "note": "Imposed starting pose for packing; not a predicted membrane orientation. Vary in replicas.",
}
json.dump(report, open(out_pdb.replace(".pdb", ".pose.json"), "w"), indent=2)
print(json.dumps(report))
