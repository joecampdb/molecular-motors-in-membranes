"""Render an MD trajectory in Blender with Cycles on the GPU, from the .npz written by export_traj.py.

  blender --background --factory-startup --python render_md.py -- \
      --npz mm1_r1.npz --view closeup|side --out DIR [--start 0] [--end 497] [--step 1]
      [--samples 64] [--res 1920 1080] [--still N]

No add-ons needed. Atoms are spheres instanced on mesh vertices, the motor's bonds are cylinders, and
a frame_change_pre handler moves every vertex each frame.

Two view-specific tricks:
  side     lipids are drawn in three periodic copies along x, so the bilayer runs past both edges of
           the frame instead of floating as an isolated patch. The motor is not copied.
  closeup  the camera tracks the motor's centre of mass, and lipid atoms between it and the camera are
           moved far off-screen for that frame. Clipping them with the camera's near plane instead
           would slice the spheres open and show their black interiors.
"""
import os
import sys

import bpy
import numpy as np
from mathutils import Vector

argv = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else []


def opt(name, default=None, n=1, cast=str):
    if name not in argv:
        return default
    i = argv.index(name)
    vals = [cast(v) for v in argv[i + 1:i + 1 + n]]
    return vals[0] if n == 1 else vals


NPZ = opt("--npz")
VIEW = opt("--view", "closeup")
OUT = opt("--out", "//render")
START = opt("--start", 0, cast=int)
END = opt("--end", -1, cast=int)
STEP = opt("--step", 1, cast=int)
SAMPLES = opt("--samples", 64, cast=int)
RES = opt("--res", [1920, 1080], n=2, cast=int)
STILL = opt("--still", None, cast=int)

d = np.load(NPZ)
POS, BOX, GRP = d["pos"], d["box"], d["grp"]
ELEM = np.array([e.decode() for e in d["elem"]])
LIG, BONDS = d["lig_idx"], d["lig_bonds"]
AXLE, NPLUS = d["axle"], int(d["nplus"])
NFRAMES = POS.shape[0]
if END < 0 or END >= NFRAMES:
    END = NFRAMES - 1

TILES = (-1, 0, 1) if VIEW == "side" else (0,)
DUMP = (0.0, 0.0, 1.0e4)        # off-camera parking spot for culled atoms
# Both views cut away whatever lipid stands between the camera and the motor. Without it the side
# view buries the subject: measured over this trajectory, 98% of the motor's atoms sit behind lipid
# in projection, because the camera looks through about 32 A of bilayer to reach it. Cutting 2 A in
# front of the motor's centre drops that to 12% while keeping 55% of the lipids, so the membrane
# still reads as a solid cross-section rather than an isolated patch.
KEEP_RADIUS = 22.0 if VIEW != "side" else None   # close-up only; < half a box length
CUT_AHEAD = 1.5 if VIEW != "side" else 2.0

# ---------------------------------------------------------------- scene
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = "CYCLES"
scene.cycles.device = "GPU"
scene.cycles.samples = SAMPLES
scene.cycles.use_denoising = True
scene.cycles.max_bounces = 8
scene.render.resolution_x, scene.render.resolution_y = RES
scene.render.resolution_percentage = 100
scene.render.image_settings.file_format = "PNG"
scene.render.image_settings.color_mode = "RGB"
scene.view_settings.view_transform = "AgX"
for look in ("AgX - Medium High Contrast", "AgX - Base Contrast", "None"):
    try:
        scene.view_settings.look = look
        break
    except TypeError:
        continue

prefs = bpy.context.preferences.addons["cycles"].preferences
for dev_type in ("OPTIX", "CUDA"):
    try:
        prefs.compute_device_type = dev_type
    except TypeError:
        continue
    prefs.refresh_devices()
    if any(dv.type == dev_type for dv in prefs.devices):
        for dv in prefs.devices:
            dv.use = (dv.type == dev_type)
        print("[setup] device %s: %s" % (dev_type, ", ".join(dv.name for dv in prefs.devices if dv.use)))
        break

# ---------------------------------------------------------------- materials
def material(name, rgba, rough=0.42, emit=0.0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    b = m.node_tree.nodes["Principled BSDF"]
    b.inputs["Base Color"].default_value = rgba
    b.inputs["Roughness"].default_value = rough
    if emit and "Emission Color" in b.inputs:
        b.inputs["Emission Color"].default_value = rgba
        b.inputs["Emission Strength"].default_value = emit
    return m


MATS = {
    # membrane: cool, desaturated, recessive
    "tail":  material("tail",  (0.205, 0.245, 0.315, 1), rough=0.68),
    "head":  material("head",  (0.330, 0.390, 0.500, 1), rough=0.55),
    "phos":  material("phos",  (0.180, 0.720, 0.880, 1), rough=0.32, emit=0.22),
    # motor: warm, saturated, faintly self-lit
    "ligC":  material("ligC",  (0.980, 0.550, 0.070, 1), rough=0.34, emit=0.10),
    "ligH":  material("ligH",  (1.000, 0.840, 0.580, 1), rough=0.42, emit=0.06),
    "ligN":  material("ligN",  (0.900, 0.350, 0.050, 1), rough=0.34, emit=0.10),
    "ligS":  material("ligS",  (1.000, 0.880, 0.250, 1), rough=0.30, emit=0.12),
    # the two atoms the measurements are about
    "axle":  material("axle",  (1.000, 0.100, 0.080, 1), rough=0.26, emit=1.10),
    "nplus": material("nplus", (0.150, 1.000, 0.450, 1), rough=0.26, emit=0.95),
    "bond":  material("bond",  (0.820, 0.440, 0.060, 1), rough=0.38),
}

# ---------------------------------------------------------------- atoms
CLOUDS = []


def cloud(name, idx, radius, mat, seg=18, ring=11, tiles=(0,), cullable=False):
    idx = np.asarray(idx, dtype=np.int64)
    n = len(idx) * len(tiles)
    me = bpy.data.meshes.new("mesh_" + name)
    me.from_pydata([(0.0, 0.0, 0.0)] * n, [], [])
    me.update()
    ob = bpy.data.objects.new("atoms_" + name, me)
    scene.collection.objects.link(ob)
    bpy.ops.mesh.primitive_uv_sphere_add(radius=radius, segments=seg, ring_count=ring)
    proto = bpy.context.object
    proto.name = "proto_" + name
    proto.data.materials.append(mat)
    for poly in proto.data.polygons:
        poly.use_smooth = True
    proto.parent = ob
    ob.instance_type = "VERTS"
    ob.show_instancer_for_render = False
    CLOUDS.append({"mesh": me, "idx": idx, "tiles": tuple(tiles), "cullable": cullable,
                   "buf": np.empty((n, 3), dtype=np.float32)})
    return CLOUDS[-1]


is_lig = np.isin(np.arange(len(GRP)), LIG)
hi = {int(AXLE[0]), int(AXLE[1]), NPLUS}
R_TAIL, R_HEAD = (0.58, 0.60) if VIEW == "side" else (0.48, 0.50)

cloud("tail", np.where((GRP == 0) & ~is_lig)[0], R_TAIL, MATS["tail"], 14, 9, TILES, True)
cloud("head", np.where((GRP == 1) & ~is_lig)[0], R_HEAD, MATS["head"], 16, 10, TILES, True)
cloud("phos", np.where(GRP == 2)[0], 1.05, MATS["phos"], 24, 14, TILES, True)
for el, rad, key in (("C", 0.62, "ligC"), ("H", 0.30, "ligH"), ("N", 0.64, "ligN"), ("S", 0.80, "ligS")):
    sel = np.array([i for i in LIG if ELEM[i] == el and int(i) not in hi], dtype=np.int64)
    if len(sel):
        cloud("lig" + el, sel, rad, MATS[key], 24, 14)
cloud("axle", np.array(sorted(int(a) for a in AXLE)), 0.95, MATS["axle"], 28, 16)
cloud("nplus", np.array([NPLUS]), 0.95, MATS["nplus"], 28, 16)
print("[setup] atoms drawn: %d (tiles %s)" % (sum(len(c["idx"]) * len(c["tiles"]) for c in CLOUDS), TILES))

# ---------------------------------------------------------------- motor bonds
bpy.ops.mesh.primitive_cylinder_add(radius=0.20, depth=1.0, vertices=14)
bond_proto = bpy.context.object
bond_proto.name = "proto_bond"
bond_proto.data.materials.append(MATS["bond"])
for poly in bond_proto.data.polygons:
    poly.use_smooth = True
bond_proto.hide_render = True
bond_proto.hide_viewport = True

BONDOBJS = []
for k, (a, b) in enumerate(BONDS):
    o = bpy.data.objects.new("bond_%03d" % k, bond_proto.data)
    o.rotation_mode = "QUATERNION"
    scene.collection.objects.link(o)
    BONDOBJS.append((o, int(a), int(b)))
print("[setup] motor bonds: %d" % len(BONDOBJS))

# ---------------------------------------------------------------- world & lights
world = bpy.data.worlds.new("world")
world.use_nodes = True
world.node_tree.nodes["Background"].inputs[0].default_value = (0.030, 0.032, 0.038, 1)
scene.world = world


def sun(name, rot, energy, color=(1.0, 1.0, 1.0), angle=0.20):
    ld = bpy.data.lights.new(name, type="SUN")
    ld.energy, ld.color, ld.angle = energy, color, angle
    o = bpy.data.objects.new(name, ld)
    o.rotation_euler = rot
    scene.collection.objects.link(o)
    return o


sun("key",  (0.92, 0.00, 0.70), 4.2)
sun("fill", (1.28, 0.00, -1.15), 1.5, (0.78, 0.84, 1.00), angle=0.6)
sun("rim",  (-0.95, 0.00, 3.05), 3.0, (1.00, 0.84, 0.68))

# ---------------------------------------------------------------- camera
cam_data = bpy.data.cameras.new("cam")
cam = bpy.data.objects.new("cam", cam_data)
cam.rotation_mode = "QUATERNION"
scene.collection.objects.link(cam)
scene.camera = cam
CAM_OFFSET = Vector((7.5, -33.0, 9.5))

if VIEW == "side":
    z_lo, z_hi = np.percentile(POS[:, :, 2], [0.3, 99.7])
    span = float(z_hi - z_lo) + 5.0
    cam_data.type = "ORTHO"
    cam_data.ortho_scale = span * RES[0] / RES[1]
    cam.location = (0.0, -180.0, float(z_lo + z_hi) * 0.5)
    cam.rotation_quaternion = Vector((0.0, 1.0, 0.0)).to_track_quat("-Z", "Y")
    print("[setup] side: ortho %.1f A wide, %.1f A tall" % (cam_data.ortho_scale, span))
else:
    cam_data.type = "PERSP"
    cam_data.lens = 50.0
    cam_data.clip_start = 0.1
    cam_data.dof.use_dof = True
    cam_data.dof.focus_distance = CAM_OFFSET.length
    cam_data.dof.aperture_fstop = 2.4
    cam.rotation_quaternion = (-CAM_OFFSET).to_track_quat("-Z", "Y")
    print("[setup] closeup: camera %.1f A from the motor" % CAM_OFFSET.length)

# unit vector pointing from the camera toward the motor
if VIEW == "side":
    VIEWDIR = np.array((0.0, 1.0, 0.0), dtype=np.float32)
else:
    VIEWDIR = np.array((-CAM_OFFSET).normalized()[:], dtype=np.float32)

# ---------------------------------------------------------------- animation
UP = Vector((0.0, 0.0, 1.0))


def apply_frame(fr):
    i = int(np.clip(fr, 0, NFRAMES - 1))
    p = POS[i]
    bx = float(BOX[i, 0])
    centre = p[LIG].mean(axis=0)

    rel = p - centre
    visible = (rel @ VIEWDIR) > -CUT_AHEAD      # drop whatever stands in front of the motor
    if KEEP_RADIUS is not None:
        visible = visible & (np.einsum("ij,ij->i", rel, rel) < KEEP_RADIUS ** 2)

    for c in CLOUDS:
        idx, buf, tiles = c["idx"], c["buf"], c["tiles"]
        m = len(idx)
        hide = ~visible[idx] if c["cullable"] else None
        any_hidden = hide is not None and bool(hide.any())
        for t, off in enumerate(tiles):
            chunk = buf[t * m:(t + 1) * m]
            chunk[:] = p[idx]
            if off:
                chunk[:, 0] += off * bx
            if any_hidden:
                chunk[hide] = DUMP
        c["mesh"].vertices.foreach_set("co", buf.reshape(-1))
        c["mesh"].update()

    for o, a, b in BONDOBJS:
        va, vb = Vector(p[a].tolist()), Vector(p[b].tolist())
        delta = vb - va
        ln = delta.length
        o.location = (va + vb) * 0.5
        o.scale = (1.0, 1.0, max(ln, 1e-4))
        o.rotation_quaternion = UP.rotation_difference(delta.normalized() if ln > 1e-6 else UP)

    if VIEW != "side":
        cam.location = Vector(centre.tolist()) + CAM_OFFSET


def on_frame(scn, _dg=None):
    apply_frame(scn.frame_current)


bpy.app.handlers.frame_change_pre.clear()
bpy.app.handlers.frame_change_pre.append(on_frame)

scene.frame_start, scene.frame_end, scene.frame_step = START, END, STEP
os.makedirs(OUT, exist_ok=True)
scene.render.use_overwrite = False
scene.render.use_placeholder = True

if STILL is not None:
    scene.frame_set(STILL)
    apply_frame(STILL)
    scene.render.filepath = os.path.join(OUT, "still_%04d" % STILL)
    bpy.ops.render.render(write_still=True)
    print("[done] still %d" % STILL)
else:
    scene.render.filepath = os.path.join(OUT, "frame_")
    scene.frame_set(START)
    apply_frame(START)
    print("[render] %s frames %d..%d step %d, %d samples, %dx%d"
          % (VIEW, START, END, STEP, SAMPLES, RES[0], RES[1]))
    bpy.ops.render.render(animation=True)
    print("[done] animation")
