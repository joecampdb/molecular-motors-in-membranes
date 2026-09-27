"""How much of the GAFF2 profile along each scan is torsion terms and how much is intramolecular
electrostatics? Both scans are in vacuum, where a charged tail folding onto the ring can produce a
barrier that has nothing to do with the torsion parameters."""
import json, numpy as np, openmm as mm, openmm.app as app, openmm.unit as u, parmed as pmd

def read_scan(path):
    lines = open(path).read().splitlines(); xyz = []; i = 0
    while i < len(lines):
        n = int(lines[i].strip())
        xyz.append([[float(v) for v in lines[i+2+k].split()[1:4]] for k in range(n)])
        i += n + 2
    return np.array(xyz)

p = pmd.load_file("/root/work/routeA/lig/MM1/lig.prmtop", "/root/work/routeA/lig/MM1/lig.inpcrd")
scans = np.load("/root/work/qm/scans.npy", allow_pickle=True)
nplus = [a.idx for a in p.atoms if a.name == "N2"][0]
ring_idx = [a.idx for a in p.atoms if a.name in ("C1","C2","C3","C4","C5","C6")]
out = {}
for row in scans:
    key, q = row[0], [int(v) for v in row[1:]]
    xyz = read_scan("/root/work/qm/%s/xtbscan.log" % key)
    sysm = p.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, rigidWater=False)
    for f in sysm.getForces():
        n = f.__class__.__name__
        f.setForceGroup(1 if n == "NonbondedForce" else (2 if n == "PeriodicTorsionForce" else 0))
    ctx = mm.Context(sysm, mm.VerletIntegrator(0.001*u.picoseconds), mm.Platform.getPlatformByName("CPU"))
    tot, nbe, tor, fold = [], [], [], []
    for c in xyz:
        ctx.setPositions(c*0.1)
        g = lambda grp: ctx.getState(getEnergy=True, groups=grp).getPotentialEnergy().value_in_unit(u.kilocalorie_per_mole)
        tot.append(g({0,1,2})); nbe.append(g({1})); tor.append(g({2}))
        fold.append(float(np.linalg.norm(c[nplus] - c[ring_idx].mean(axis=0))))
    r = lambda v: float(np.max(v) - np.min(v))
    out[key] = {"range_total_kcal": r(tot), "range_nonbonded_kcal": r(nbe), "range_torsion_kcal": r(tor),
                "ammonium_to_ring_A_min": round(min(fold),2), "ammonium_to_ring_A_max": round(max(fold),2)}
    print("%-6s single points on QM geometries: total span %6.1f | nonbonded %6.1f | torsion terms %5.1f"
          % (key, out[key]["range_total_kcal"], out[key]["range_nonbonded_kcal"], out[key]["range_torsion_kcal"]))
    print("       ammonium-to-ring distance varies %.1f - %.1f A along the scan" % (min(fold), max(fold)))
json.dump(out, open("/root/work/qm/decomposition.json","w"), indent=1)
