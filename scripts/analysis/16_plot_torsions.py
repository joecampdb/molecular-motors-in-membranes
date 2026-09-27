"""Overlay the GFN2-xTB and GAFF2 torsion profiles.

Usage: python 16_plot_torsions.py <torsion-validation.json> <decomposition.json> <out_prefix>
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

val = json.load(open(sys.argv[1]))
dec = json.load(open(sys.argv[2]))
out = sys.argv[3]

QM_C, FF_C = "#2a78d6", "#eb6834"
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 10,
                     "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.labelcolor": INK2, "text.color": INK,
                     "axes.titlecolor": INK, "axes.spines.top": False, "axes.spines.right": False,
                     "figure.facecolor": SURF, "axes.facecolor": SURF})

TITLES = {"twist": "Twist of the overcrowded C=C, the motor's own axle",
          "arylN": "Aryl–amine bond that orients the anchoring tail"}
NOTES = {"twist": "torsion terms dominate the force-field profile, so this is a fair test",
         "arylN": "force-field profile is dominated by vacuum electrostatics, not torsion terms"}

keys = [k for k in ("twist", "arylN") if k in val["torsions"]]
fig, axes = plt.subplots(1, len(keys), figsize=(11.6, 4.5), dpi=100)
axes = np.atleast_1d(axes)
fig.subplots_adjust(left=0.07, right=0.985, top=0.80, bottom=0.14, wspace=0.22)

for ax, k in zip(axes, keys):
    t = val["torsions"][k]
    x = np.array(t["angle_deg"])
    ax.plot(x, t["qm_kcal"], color=QM_C, lw=2.0, solid_capstyle="round", label="GFN2-xTB (reference)")
    ax.plot(x, t["ff_relaxed_kcal"], color=FF_C, lw=2.0, solid_capstyle="round", label="GAFF2 (used in the MD)")
    for s, c in ((t["stiffness_qm"], QM_C), (t["stiffness_ff"], FF_C)):
        if s:
            ax.axvline(s["minimum_deg"], color=c, lw=0.9, alpha=0.35)
    sq, sf = t["stiffness_qm"], t["stiffness_ff"]
    if sq and sf and k == "twist":
        # only annotate the fair comparison, and use the shortest arc between the two minima
        off = (sf["minimum_deg"] - sq["minimum_deg"] + 180.0) % 360.0 - 180.0
        y = ax.get_ylim()[1] * 0.06
        ax.annotate("", xy=(sf["minimum_deg"], y), xytext=(sq["minimum_deg"], y),
                    arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.0))
        ax.text((sq["minimum_deg"] + sf["minimum_deg"]) / 2, y * 1.5,
                "%.0f° offset in the resting twist" % abs(off),
                ha="center", va="bottom", color=INK2, fontsize=9)
    d = dec.get(k, {})
    ax.set_title(TITLES[k], loc="left", fontsize=11.5, pad=16)
    ax.text(0, 1.02, NOTES[k], transform=ax.transAxes, fontsize=9, color=MUTED, va="bottom")
    ax.set_xlabel("dihedral (degrees)")
    ax.grid(axis="y", color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    ax.set_xlim(x.min(), x.max())
axes[0].set_ylabel("relative energy (kcal/mol)")
axes[0].legend(frameon=False, fontsize=9.5, loc="upper center", handlelength=1.8)

fig.text(0.07, 0.945, "Does the force field describe the motor's own core?", fontsize=13,
         color=INK, weight="semibold")
fig.text(0.07, 0.905,
         "Relaxed scans of the same molecule, same atom indexing, both in vacuum. Curvature is right; the resting geometry is not.",
         fontsize=9.5, color=INK2)
for ext in ("png", "svg"):
    fig.savefig("%s.%s" % (out, ext), dpi=200 if ext == "png" else None, facecolor=SURF)
print("wrote %s.png / .svg" % out)
