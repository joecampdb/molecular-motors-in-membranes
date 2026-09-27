"""Convert the packmol-memgen/tleap Amber system to GROMACS with ParmEd and audit it.

Usage: python 04_to_gromacs.py <sys.prmtop> <sys.inpcrd|rst7> <outprefix>
Writes <outprefix>.top, <outprefix>.gro and <outprefix>.audit.json.

Why the OpenMM round trip: this prmtop carries non-uniform 1-4 Lennard-Jones scaling (SCNB 2.0 and
6.0). ParmEd 4.3.1's direct Amber->GROMACS writer then emits "[ defaults ] ... gen-pairs yes fudgeLJ 1"
with bare 3-column [ pairs ], which makes GROMACS apply UNscaled 1-4 LJ (measured: LJ-14 doubled,
+4.9 MJ/mol per system). Building the ParmEd structure from the OpenMM System instead gives every 1-4
pair explicit parameters, so GROMACS reproduces Amber's per-pair scaling exactly.
"""
import collections
import json
import re
import sys

import openmm.app as app
import openmm.unit as u
import parmed as pmd

prmtop, inpcrd, out = sys.argv[1:4]
p = pmd.load_file(prmtop, inpcrd)

scee = sorted(set(round(t.scee, 4) for t in p.dihedral_types))
scnb = sorted(set(round(t.scnb, 4) for t in p.dihedral_types))
nonstandard = collections.Counter()
for d in p.dihedrals:
    if d.type is not None and not d.ignore_end and round(d.type.scnb, 3) != 2.0:
        nonstandard[(round(d.type.scnb, 3), d.atom1.residue.name, d.atom1.type, d.atom2.type, d.atom3.type, d.atom4.type)] += 1

system = p.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=1.0 * u.nanometer, constraints=None, rigidWater=False)
s = pmd.openmm.load_topology(p.topology, system, xyz=p.positions)
s.box = p.box
s.save(out + ".top", format="gromacs", overwrite=True)
s.save(out + ".gro", overwrite=True)

# Rename ParmEd's generic multi-residue molecule names (systemN) to POPE/POPG by head-group residue.
txt = open(out + ".top").read()
blocks = re.split(r"(?=\[ moleculetype \])", txt)
rename = {}
for b in blocks:
    m = re.match(r"\[ moleculetype \]\s*\n(?:;[^\n]*\n)*\s*(\S+)", b)
    if not m or not m.group(1).startswith("system"):
        continue
    resnames = set(re.findall(r"^\s*\d+\s+\S+\s+\d+\s+(\S+)\s+\S+\s+\d+", b, flags=re.M))
    if "PE" in resnames and "PGR" not in resnames:
        rename[m.group(1)] = "POPE"
    elif "PGR" in resnames and "PE" not in resnames:
        rename[m.group(1)] = "POPG"
for old, new in rename.items():
    txt = re.sub(rf"\b{old}\b", new, txt)
open(out + ".top", "w").write(txt)

# Audit the GROMACS-facing facts textually, independent of ParmEd's reader.
defaults = re.search(r"\[ defaults \]\s*\n(?:;[^\n]*\n)*\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", txt)
pair_cols = collections.Counter()
in_pairs = False
for line in txt.splitlines():
    if line.startswith("["):
        in_pairs = line.strip() == "[ pairs ]"
        continue
    if in_pairs and line.strip() and not line.startswith(";"):
        pair_cols[len(line.split())] += 1
res_counts = collections.Counter(r.name for r in p.residues)
audit = {
    "prmtop": prmtop, "coordinates": inpcrd, "atoms": len(p.atoms), "residues": len(p.residues),
    "net_charge_e": round(float(sum(a.charge for a in p.atoms)), 6),
    "ligand_charge_e": round(float(sum(a.charge for r in p.residues if r.name == "LIG" for a in r.atoms)), 6),
    "waters": res_counts.get("WAT", 0), "ions": {k: v for k, v in res_counts.items() if k in ("Na+", "Cl-", "K+")},
    "POPE_total": res_counts.get("PE", 0), "POPG_total": res_counts.get("PGR", 0),
    "lipids_total": res_counts.get("PE", 0) + res_counts.get("PGR", 0),
    "box_A": [round(float(x), 3) for x in p.box[:3]],
    "amber_scee_values": scee, "amber_scnb_values": scnb,
    "dihedrals_with_nonstandard_scnb": [{"scnb": k[0], "residue": k[1], "types": list(k[2:]), "count": v} for k, v in sorted(nonstandard.items())],
    "gromacs_defaults": {"nbfunc": defaults.group(1), "comb_rule": defaults.group(2), "gen_pairs": defaults.group(3), "fudgeLJ": defaults.group(4), "fudgeQQ": defaults.group(5)} if defaults else None,
    "gromacs_pairs_lines_by_column_count": dict(pair_cols),
    "gromacs_pairs_all_explicit": all(c >= 5 for c in pair_cols) if pair_cols else False,
    "gromacs_top_has_settles": "[ settles ]" in txt,
    "gromacs_top_moleculetypes": txt.count("[ moleculetype ]"),
    "moleculetype_renames": rename,
    "residue_counts": dict(res_counts),
}
json.dump(audit, open(out + ".audit.json", "w"), indent=2)
print(json.dumps({k: v for k, v in audit.items() if k not in ("residue_counts",)}))
