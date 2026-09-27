"""Build a render-ready .npz for Blender from the GROMACS render trajectory.

  python export_traj.py <render.gro> <render.xtc> <lig.mol2> <out.npz> [smooth_frames]

Everything expensive or subtle happens here, so the Blender script only has to draw:

1. Selection. The motor keeps every atom; lipids keep heavy atoms only (their hydrogens are 59% of
   all atoms and read as pure speckle). Water and ions are dropped.
2. Time-unwrapping, per residue. `trjconv -pbc mol -center` re-wraps whole molecules as the centring
   target moves, so 4.3% of atoms jump more than 10 A between neighbouring frames, up to a full box
   length. On screen those atoms teleport. Unwrapping each residue against the previous frame with
   the minimum image removes it.
3. Smoothing. A centred running mean over `smooth_frames` frames (200 ps each). At 5 frames this is
   a 1 ns window and removes about 70% of the frame-to-frame jitter, which is what makes the motion
   readable. It is a display choice and must be stated wherever the video is shown.
4. Centring. The motor's centre of mass goes to x = y = 0 every frame and the mean phosphorus plane
   to z = 0, so the motor holds one screen position and the cameras can be placed absolutely.
5. Re-wrapping, per lipid residue, to its nearest image of the motor in x and y. This keeps the patch
   compact around the subject. A residue only ever flips at half a box length, which the side view's
   periodic tiling hides and the close-up's 22 A cull radius stays inside.

Arrays: pos (nframes, natoms, 3) float32 Angstrom, box (nframes, 3), elem (natoms) 'S2',
grp (natoms) uint8 0=tail 1=head 2=phosphorus 3=motor, lig_idx, lig_bonds, axle (2,), nplus ().
"""
import sys

import numpy as np
import MDAnalysis as mda

gro, xtc, mol2, out = sys.argv[1:5]
SMOOTH = int(sys.argv[5]) if len(sys.argv) > 5 else 5

u = mda.Universe(gro, xtc)
LIPIDS = ("PA", "OL", "PE", "PGR")

names_all = np.array([a.name for a in u.atoms])
resn_all = np.array([a.resname for a in u.atoms])
is_h = np.array([n[0] == "H" for n in names_all])
keep = np.where((resn_all == "LIG") | (np.isin(resn_all, LIPIDS) & ~is_h))[0]

sub = u.atoms[keep]
names = names_all[keep]
resn = resn_all[keep]
resix = sub.resindices.copy()
n = len(keep)

elem = np.array([nm[0] for nm in names], dtype="S2")
grp = np.zeros(n, dtype=np.uint8)
grp[np.isin(resn, ["PE", "PGR"])] = 1
grp[names == "P31"] = 2
lig_mask = resn == "LIG"
grp[lig_mask] = 3
lig_idx = np.where(lig_mask)[0].astype(np.int32)
assert len(lig_idx) == 65, len(lig_idx)
dropped = len(u.atoms) - n
print("kept %d of %d atoms (dropped %d: lipid H, water, ions)" % (n, len(u.atoms), dropped))
print("  tails %d, heads %d, phosphorus %d, motor %d"
      % ((grp == 0).sum(), (grp == 1).sum(), (grp == 2).sum(), (grp == 3).sum()))

bonds = []
sec = None
for line in open(mol2):
    if line.startswith("@<TRIPOS>"):
        sec = line.strip()
        continue
    if sec == "@<TRIPOS>BOND":
        p = line.split()
        if len(p) >= 4:
            bonds.append((int(p[1]) - 1, int(p[2]) - 1))
lig_bonds = np.array([[lig_idx[a], lig_idx[b]] for a, b in bonds], dtype=np.int32)

nf = len(u.trajectory)
pos = np.empty((nf, n, 3), dtype=np.float32)
box = np.empty((nf, 3), dtype=np.float32)
for i, ts in enumerate(u.trajectory):
    pos[i] = ts.positions[keep]
    box[i] = ts.dimensions[:3]

# --- 2. unwrap each residue in time -------------------------------------------------------------
order = np.argsort(resix, kind="stable")
_, starts = np.unique(resix[order], return_index=True)
res_groups = np.split(order, starts[1:])
raw_jump = float((np.linalg.norm(np.diff(pos[:60], axis=0), axis=2) > 10).mean() * 100)
for i in range(1, nf):
    L = box[i]
    for g in res_groups:
        delta = pos[i, g].mean(axis=0) - pos[i - 1, g].mean(axis=0)
        shift = np.round(delta / L) * L
        if shift.any():
            pos[i, g] -= shift
post_jump = float((np.linalg.norm(np.diff(pos[:60], axis=0), axis=2) > 10).mean() * 100)
print("atoms jumping >10 A per frame: %.2f%% raw -> %.3f%% unwrapped" % (raw_jump, post_jump))

# --- 3. centred running mean --------------------------------------------------------------------
if SMOOTH > 1:
    pad = SMOOTH // 2
    padded = np.concatenate([pos[:1]] * pad + [pos] + [pos[-1:]] * pad, axis=0)
    csum = np.cumsum(padded, axis=0, dtype=np.float64)
    csum = np.concatenate([np.zeros((1, n, 3)), csum], axis=0)
    pos = ((csum[SMOOTH:] - csum[:-SMOOTH]) / SMOOTH).astype(np.float32)[:nf]
    print("smoothed over %d frames (%.1f ns window); median step now %.2f A"
          % (SMOOTH, SMOOTH * 0.2, np.median(np.linalg.norm(np.diff(pos[:60], axis=0), axis=2))))

# --- 4. centre on the motor, and on the phosphorus plane in z ------------------------------------
lig_com = pos[:, lig_idx].mean(axis=1)
pos[:, :, 0] -= lig_com[:, 0][:, None]
pos[:, :, 1] -= lig_com[:, 1][:, None]
pos[:, :, 2] -= float(pos[:, grp == 2, 2].mean())

# --- 5. nearest image of each lipid residue to the motor, in x and y only ------------------------
lip_groups = [g for g in res_groups if grp[g[0]] != 3]
for i in range(nf):
    L = box[i]
    for g in lip_groups:
        c = pos[i, g].mean(axis=0)
        for ax in (0, 1):
            k = np.round(c[ax] / L[ax])
            if k:
                pos[i, g, ax] -= k * L[ax]

np.savez_compressed(out, pos=pos, box=box, elem=elem, grp=grp, lig_idx=lig_idx,
                    lig_bonds=lig_bonds, axle=lig_idx[[6, 14]].astype(np.int32),
                    nplus=np.int32(lig_idx[31]), smooth=np.int32(SMOOTH))
print("frames %d  atoms %d  bonds %d  -> %s" % (nf, n, len(lig_bonds), out))
print("axle atoms %s %s, ammonium %s" % (names[lig_idx[6]], names[lig_idx[14]], names[lig_idx[31]]))
