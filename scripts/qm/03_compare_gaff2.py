"""Compare the GAFF2 torsion profile against a GFN2-xTB relaxed scan, for the same molecule and the
same atom indexing.

Usage: python ff_scan.py <qm_dir> <out.json>

For every point of the QM scan: restrain the same dihedral at the same value, minimise under GAFF2,
and record the force-field energy with the restraint excluded. A single-point GAFF2 energy on the
unrelaxed QM geometry is also recorded, which separates "the force field prefers a different
geometry" from "the force field prefers a different energy".

Both curves are then referenced to their own minimum, so what is compared is the shape.

The stiffness near the minimum is turned into the quantity that actually matters for a simulation:
the root-mean-square thermal fluctuation of that dihedral at 300 K, sqrt(kT/k), from a quadratic fit
over the innermost points.
"""
import glob
import json
import os
import sys

import numpy as np
import openmm as mm
import openmm.app as app
import openmm.unit as u
import parmed as pmd

QM, OUT = sys.argv[1], sys.argv[2]
HARTREE = 627.5094740631      # kcal/mol
RT = 0.001987204259 * 300.0   # kcal/mol at 300 K


def read_scan(path):
    """Read an xtb multi-frame xyz; return (energies in Hartree, coords in Angstrom)."""
    lines = open(path).read().splitlines()
    e, xyz, i = [], [], 0
    while i < len(lines):
        n = int(lines[i].strip())
        tok = [t for t in lines[i + 1].replace(":", " ").split() if _isfloat(t)]
        e.append(float(tok[0]))
        xyz.append([[float(v) for v in lines[i + 2 + k].split()[1:4]] for k in range(n)])
        i += n + 2
    return np.array(e), np.array(xyz)


def _isfloat(t):
    try:
        float(t)
        return True
    except ValueError:
        return False


def dihedral(c, q):
    a, b, cc, d = [c[i] for i in q]
    b0, b1, b2 = a - b, cc - b, d - cc
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def stiffness(x_deg, y_kcal, half_width=25.0):
    """Quadratic fit about the minimum -> force constant and RMS fluctuation at 300 K."""
    x = np.asarray(x_deg, float)
    y = np.asarray(y_kcal, float)
    x0 = x[int(np.argmin(y))]
    m = np.abs(x - x0) <= half_width
    if m.sum() < 4:
        return None
    a = np.polyfit(np.radians(x[m] - x0), y[m], 2)[0]      # y = a * dtheta^2 + ...
    k = 2.0 * a                                            # kcal/mol/rad^2
    if k <= 0:
        return None
    return {"minimum_deg": float(x0), "k_kcal_per_mol_rad2": float(k),
            "rms_fluctuation_deg_300K": float(np.degrees(np.sqrt(RT / k)))}


p = pmd.load_file("/root/work/routeA/lig/MM1/lig.prmtop", "/root/work/routeA/lig/MM1/lig.inpcrd")
scans = np.load(os.path.join(QM, "scans.npy"), allow_pickle=True)
result = {"method_qm": "GFN2-xTB relaxed scan, charge +1",
          "method_ff": "GAFF2 / AM1-BCC, vacuum, relaxed under the same dihedral restraint",
          "note_qm": ("GFN2-xTB is semi-empirical tight binding. The profile near the minimum is the "
                      "trustworthy part; a strongly twisted overcrowded alkene acquires diradical "
                      "character that this method does not describe well."),
          "torsions": {}}

for row in scans:
    key, q = row[0], [int(v) for v in row[1:]]
    e_qm, xyz = read_scan(os.path.join(QM, key, "xtbscan.log"))
    angles = np.array([dihedral(c, q) for c in xyz])

    system = p.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    for f in system.getForces():
        f.setForceGroup(0)
    rest = mm.CustomTorsionForce("0.5*k*dmin^2; dmin=min(d, 6.2831853-d); d=abs(theta-theta0)")
    rest.addPerTorsionParameter("k")
    rest.addPerTorsionParameter("theta0")
    rest.addTorsion(q[0], q[1], q[2], q[3], [2000.0, 0.0])
    rest.setForceGroup(1)
    system.addForce(rest)
    plat = mm.Platform.getPlatformByName("CPU")
    ctx = mm.Context(system, mm.VerletIntegrator(0.001 * u.picoseconds), plat)

    e_ff_relaxed, e_ff_single, ang_ff = [], [], []
    for c, target in zip(xyz, angles):
        ctx.setPositions(c * 0.1)                                    # A -> nm
        rest.setTorsionParameters(0, q[0], q[1], q[2], q[3], [0.0, np.radians(target)])
        rest.updateParametersInContext(ctx)
        e_ff_single.append(ctx.getState(getEnergy=True, groups={0}).getPotentialEnergy()
                           .value_in_unit(u.kilocalorie_per_mole))
        rest.setTorsionParameters(0, q[0], q[1], q[2], q[3], [2000.0, np.radians(target)])
        rest.updateParametersInContext(ctx)
        mm.LocalEnergyMinimizer.minimize(ctx, tolerance=1e-4 * u.kilojoule_per_mole / u.nanometer,
                                         maxIterations=2000)
        st = ctx.getState(getEnergy=True, getPositions=True, groups={0})
        e_ff_relaxed.append(st.getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole))
        ang_ff.append(dihedral(np.array(st.getPositions().value_in_unit(u.angstrom)), q))

    qm = (e_qm - e_qm.min()) * HARTREE
    ffr = np.array(e_ff_relaxed) - np.min(e_ff_relaxed)
    ffs = np.array(e_ff_single) - np.min(e_ff_single)
    order = np.argsort(angles)
    inner = np.abs(angles - angles[int(np.argmin(qm))]) <= 40.0
    result["torsions"][key] = {
        "atoms_0based": q, "n_points": len(angles),
        "angle_deg": angles[order].round(2).tolist(),
        "qm_kcal": qm[order].round(3).tolist(),
        "ff_relaxed_kcal": ffr[order].round(3).tolist(),
        "ff_single_point_on_qm_geometry_kcal": ffs[order].round(3).tolist(),
        "qm_range_kcal": float(qm.max()),
        "ff_range_kcal": float(ffr.max()),
        "rmsd_ff_vs_qm_kcal": float(np.sqrt(np.mean((ffr - qm) ** 2))),
        "rmsd_ff_vs_qm_inner40deg_kcal": float(np.sqrt(np.mean((ffr[inner] - qm[inner]) ** 2))),
        "max_abs_error_kcal": float(np.max(np.abs(ffr - qm))),
        "stiffness_qm": stiffness(angles, qm),
        "stiffness_ff": stiffness(angles, ffr),
    }
    t = result["torsions"][key]
    print("\n%s  (atoms %s)" % (key, q))
    print("  QM span %.1f kcal/mol over %+.0f..%+.0f deg; GAFF2 span %.1f"
          % (t["qm_range_kcal"], min(angles), max(angles), t["ff_range_kcal"]))
    print("  GAFF2 vs QM: RMSD %.2f kcal/mol overall, %.2f within 40 deg of the minimum, worst %.2f"
          % (t["rmsd_ff_vs_qm_kcal"], t["rmsd_ff_vs_qm_inner40deg_kcal"], t["max_abs_error_kcal"]))
    for lab, s in (("QM   ", t["stiffness_qm"]), ("GAFF2", t["stiffness_ff"])):
        if s:
            print("  %s minimum %+6.1f deg, k = %7.1f kcal/mol/rad^2, RMS fluctuation at 300 K = %4.1f deg"
                  % (lab, s["minimum_deg"], s["k_kcal_per_mol_rad2"], s["rms_fluctuation_deg_300K"]))

json.dump(result, open(OUT, "w"), indent=1)
print("\nwrote", OUT)
