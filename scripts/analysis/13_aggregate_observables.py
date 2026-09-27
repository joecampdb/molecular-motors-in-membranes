"""Aggregate per-replica observables (from 12_observables_all_gmx.sh) into per-motor statistics.

Usage: python 13_aggregate_observables.py <work_root> <out_dir> [lipids_per_leaflet=44.5]
Reads <work_root>/sys/<M>/prod/r<N>/obs/{series.dat,box.xvg}; copies the series into <out_dir>/series/
and writes <out_dir>/summary.json and <out_dir>/summary.md.
Statistics: per replica, means over the analysis window (30-100 ns primary; 10-100 ns for comparison);
per motor, mean and standard error of the mean across the five replica means (n = 5), which is the
only honest uncertainty here since frames within a replica are strongly autocorrelated.
"""
import json
import os
import shutil
import sys

import numpy as np

root, out = sys.argv[1], sys.argv[2]
nleaf = float(sys.argv[3]) if len(sys.argv) > 3 else 44.5
os.makedirs(os.path.join(out, "series"), exist_ok=True)
WINDOWS = {"30-100 ns": (30000, 100001), "10-100 ns": (10000, 100001)}
COLS = ["tilt_deg", "axle_z_A", "nplus_z_A", "pp_nm"]
summary = {"lipids_per_leaflet": nleaf, "frames_per_replica": None, "windows": list(WINDOWS), "motors": {}}
md = []
for M in ("MM1", "MM2"):
    reps = {}
    for R in range(1, 6):
        d = os.path.join(root, "sys", M, "prod", f"r{R}", "obs")
        s = np.loadtxt(os.path.join(d, "series.dat"))
        box = np.loadtxt(os.path.join(d, "box.xvg"), comments=("#", "@"))
        apl = box[:, 1] * box[:, 2] * 100.0 / nleaf  # nm^2 -> A^2 per lipid
        shutil.copy(os.path.join(d, "series.dat"), os.path.join(out, "series", f"{M}_r{R}_series.dat"))
        shutil.copy(os.path.join(d, "box.xvg"), os.path.join(out, "series", f"{M}_r{R}_box.xvg"))
        summary["frames_per_replica"] = int(s.shape[0])
        rep = {"windows": {}, "blocks_10ns": [], "in_membrane": {}}
        t = s[:, 0]
        for wname, (a, b) in WINDOWS.items():
            m = (t >= a) & (t < b)
            mb = (box[:, 0] >= a) & (box[:, 0] < b)
            w = {c: float(s[m, i + 1].mean()) for i, c in enumerate(COLS)}
            w["tilt_sd_deg"] = float(s[m, 1].std())
            w["frac_tilt_lt_30"] = float((s[m, 1] < 30).mean())
            w["apl_A2"] = float(apl[mb].mean())
            rep["windows"][wname] = w
        for b0 in range(0, 100000, 10000):
            m = (t >= b0) & (t < b0 + 10000)
            rep["blocks_10ns"].append({"start_ns": b0 / 1000, "tilt_deg": float(s[m, 1].mean()), "axle_z_A": float(s[m, 2].mean()), "pp_nm": float(s[m, 4].mean())})
        rep["in_membrane"] = {"axle_z_min_A": float(s[:, 2].min()), "axle_z_max_A": float(s[:, 2].max()),
                              "nplus_z_min_A": float(s[:, 3].min()), "nplus_z_max_A": float(s[:, 3].max())}
        rep["tilt_hist_30_100_ns"] = np.histogram(s[(t >= 30000), 1], bins=np.arange(0, 91, 10))[0].tolist()
        reps[f"r{R}"] = rep
    motor = {"replicas": reps, "across_replicas": {}}
    for wname in WINDOWS:
        acc = {}
        for key in COLS + ["tilt_sd_deg", "frac_tilt_lt_30", "apl_A2"]:
            vals = np.array([reps[r]["windows"][wname][key] for r in reps])
            acc[key] = {"mean": float(vals.mean()), "sem": float(vals.std(ddof=1) / np.sqrt(len(vals))), "replica_values": [round(float(v), 3) for v in vals]}
        motor["across_replicas"][wname] = acc
    pooled = np.sum([reps[r]["tilt_hist_30_100_ns"] for r in reps], axis=0)
    motor["pooled_tilt_hist_30_100_ns"] = {"bins_deg": list(range(0, 91, 10)), "fraction": (pooled / pooled.sum()).round(4).tolist()}
    summary["motors"][M] = motor
    # markdown
    md.append(f"\n### {M}\n")
    md.append("| Replica | Tilt 30-100 ns (deg, mean +- sd) | Axle height above P midplane (A) | N+ height (A) | P-P thickness (nm) | APL (A^2) | fraction tilt < 30 deg |")
    md.append("| --- | --- | --- | --- | --- | --- | --- |")
    for r, rep in reps.items():
        w = rep["windows"]["30-100 ns"]
        md.append(f"| {r} | {w['tilt_deg']:.1f} +- {w['tilt_sd_deg']:.1f} | {w['axle_z_A']:.1f} | {w['nplus_z_A']:.1f} | {w['pp_nm']:.2f} | {w['apl_A2']:.1f} | {w['frac_tilt_lt_30']:.2f} |")
    a = motor["across_replicas"]["30-100 ns"]
    md.append(f"| **mean +- SEM (n=5)** | **{a['tilt_deg']['mean']:.1f} +- {a['tilt_deg']['sem']:.1f}** | **{a['axle_z_A']['mean']:.1f} +- {a['axle_z_A']['sem']:.1f}** | **{a['nplus_z_A']['mean']:.1f} +- {a['nplus_z_A']['sem']:.1f}** | **{a['pp_nm']['mean']:.2f} +- {a['pp_nm']['sem']:.2f}** | **{a['apl_A2']['mean']:.1f} +- {a['apl_A2']['sem']:.1f}** | **{a['frac_tilt_lt_30']['mean']:.2f}** |")
    h = motor["pooled_tilt_hist_30_100_ns"]["fraction"]
    md.append("\nPooled tilt histogram, 30-100 ns, 10 deg bins from 0: " + ", ".join(f"{100*f:.0f}%" for f in h))
json.dump(summary, open(os.path.join(out, "summary.json"), "w"), indent=1)
open(os.path.join(out, "summary_tables.md"), "w").write("\n".join(md) + "\n")
print("\n".join(md))
