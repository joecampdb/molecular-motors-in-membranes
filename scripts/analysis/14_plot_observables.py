"""Summary figure for the ten production replicas: tilt time courses per motor, per-replica means
against the paper's approximate comparators, and axle height against the phosphate plane.

Usage: python 14_plot_observables.py <analysis_dir>
Reads <analysis_dir>/series/<M>_r<N>_series.dat and summary.json; writes tilt_depth_summary.{svg,png}.
Colour: replicas 1-5 use categorical slots 1-5 of the validated reference palette in fixed order;
motor means and references are drawn in ink, never in a series colour.
"""
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

A = sys.argv[1]
S = json.load(open(f"{A}/summary.json"))
SER = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]   # slots 1-5, light mode
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
plt.rcParams.update({"font.family": ["Segoe UI", "DejaVu Sans", "sans-serif"], "font.size": 9.5,
                     "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "xtick.color": MUTED, "ytick.color": MUTED,
                     "axes.labelcolor": INK2, "text.color": INK, "axes.titlecolor": INK,
                     "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": SURF, "axes.facecolor": SURF})
fig = plt.figure(figsize=(11, 7.2), dpi=100)
gs = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.5, wspace=0.28, left=0.07, right=0.98, top=0.84, bottom=0.08)
PAPER = {"MM1": 15, "MM2": 60}


def running_mean(y, w=100):
    c = np.cumsum(np.insert(y, 0, 0.0))
    return (c[w:] - c[:-w]) / w


# Row 1: tilt vs time, one panel per motor, replicas in fixed hues
axes_t = []
for j, M in enumerate(("MM1", "MM2")):
    ax = fig.add_subplot(gs[0, j]); axes_t.append(ax)
    ends = []
    for r in range(1, 6):
        s = np.loadtxt(f"{A}/series/{M}_r{r}_series.dat")
        y = running_mean(s[:, 1]); t = s[99:, 0] / 1000.0
        ax.plot(t, y, color=SER[r - 1], lw=1.5, solid_capstyle="round", label=f"replica {r}")
        ends.append([y[-1], r])
    # selective direct labels at the line ends, spread apart if they collide
    ends.sort()
    for k in range(1, len(ends)):
        if ends[k][0] - ends[k - 1][0] < 4.5:
            ends[k][0] = ends[k - 1][0] + 4.5
    for yv, r in ends:
        ax.text(101.2, yv, f"r{r}", color=INK2, fontsize=8.5, va="center")
    ax.axvspan(0, 30, color=GRID, alpha=0.45, lw=0)
    ax.text(15, 86, "excluded\n(settling)", ha="center", va="top", color=MUTED, fontsize=8)
    ax.axhline(PAPER[M], color=INK2, lw=0.8, alpha=0.35)
    ax.text(31.5, PAPER[M] + 1.2, f"paper approx. {PAPER[M]}°", ha="left", va="bottom", color=MUTED, fontsize=8,
            bbox=dict(boxstyle="round,pad=0.15", fc=SURF, ec="none", alpha=0.9))
    ax.set_xlim(0, 108); ax.set_ylim(0, 90); ax.set_yticks(range(0, 91, 15))
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
    ax.set_title(f"{M}: axle tilt from the membrane plane, 1 ns running mean", loc="left", fontsize=10.5, pad=8)
    ax.set_xlabel("time (ns)")
    if j == 0:
        ax.set_ylabel("tilt from plane (degrees)")
h, l = axes_t[0].get_legend_handles_labels()
fig.legend(h, l, frameon=False, ncol=5, loc="upper left", bbox_to_anchor=(0.07, 0.905), fontsize=8.5, handlelength=1.6, columnspacing=1.2)

# Row 2: per-replica means (30-100 ns) with motor mean +- SEM, tilt (left) and axle height (right)
def dot_panel(ax, key, ylabel, refs, ylim):
    for j, M in enumerate(("MM1", "MM2")):
        acc = S["motors"][M]["across_replicas"]["30-100 ns"][key]
        vals = acc["replica_values"]
        for r, v in enumerate(vals, start=1):
            ax.scatter(j + (r - 3) * 0.07, v, s=64, color=SER[r - 1], edgecolor=SURF, linewidth=1.5, zorder=3)
        m, e = acc["mean"], acc["sem"]
        ax.errorbar(j + 0.36, m, yerr=e, fmt="_", color=INK, ms=14, mew=1.6, capsize=4, elinewidth=1.2, zorder=4)
        ax.text(j + 0.44, m, f"{m:.1f} ± {e:.1f}", va="center", ha="left", color=INK, fontsize=9)
    for label, (y, kind) in refs.items():
        if kind == "line":
            ax.axhline(y, color=INK2, lw=0.8, alpha=0.35); ax.text(1.72, y + 0.6, label, ha="right", va="bottom", color=MUTED, fontsize=8)
        else:
            for j, M in enumerate(("MM1", "MM2")):
                ax.scatter(j - 0.36, y[M], marker="D", s=46, facecolor=SURF, edgecolor=INK2, linewidth=1.2, zorder=3)
            ax.text(-0.36, y["MM2"] + 3.5, label, ha="center", va="bottom", color=MUTED, fontsize=8)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["MM1", "MM2"], color=INK); ax.set_xlim(-0.6, 1.75); ax.set_ylim(*ylim)
    ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.set_ylabel(ylabel)

ax3 = fig.add_subplot(gs[1, 0])
dot_panel(ax3, "tilt_deg", "mean tilt from plane, 30–100 ns (degrees)", {"paper approx. (CHARMM36)": (PAPER, "marker")}, (0, 70))
ax3.set_title("Replica means (dots) and motor mean ± SEM, n = 5 (bar)", loc="left", fontsize=10.5, pad=8)
pp = np.mean([S["motors"][M]["across_replicas"]["30-100 ns"]["pp_nm"]["mean"] for M in ("MM1", "MM2")]) * 10 / 2
ax4 = fig.add_subplot(gs[1, 1])
dot_panel(ax4, "axle_z_A", "axle midpoint height above bilayer midplane (Å)", {f"upper phosphate plane ({pp:.1f} Å)": (pp, "line")}, (8, 24))
ax4.set_title("Axle depth: both motors sit just below the phosphates", loc="left", fontsize=10.5, pad=8)
fig.text(0.07, 0.965, "MM1 and MM2 in POPE:POPG 1:3, Amber Lipid21 + GAFF2, 5 × 100 ns each", fontsize=12, color=INK, weight="semibold")
fig.text(0.07, 0.935, "Independent replication with a different force field than the source paper; replicas share one starting pose and differ by velocity seed.", fontsize=8.5, color=INK2)
for ext in ("svg", "png"):
    fig.savefig(f"{A}/tilt_depth_summary.{ext}", dpi=200 if ext == "png" else None, facecolor=SURF)
print("written", f"{A}/tilt_depth_summary.svg/png")
