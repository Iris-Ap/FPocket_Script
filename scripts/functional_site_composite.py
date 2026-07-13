"""
functional_site_composite.py -- Re-score each TF at its FUNCTIONAL-site pocket (the pocket
nearest the known drug residue) instead of the highest-druggability decoy, using the same
control baseline and weights as zscore_composite.py. Reveals the "inversion": geometry
rates the decoy druggable and the real drug site undruggable.

Functional residues: FOXA1 Cys258 (Wing2), FOXM1 His287 (DBD), FOXO3 Arg249 (DBD).
Run: python3 functional_site_composite.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import grid_enclosure as GE
import zscore_composite as Z   # reuses control baseline (Z.stats), weights, orientation

FUNCTIONAL = {"FOXA1_7VOX": ("CYS", 258), "FOXM1_3G73": ("HIS", 287), "FOXO3_2UZK": ("ARG", 249)}


def residue_atoms(name, resseq):
    p = f"{C.BASE}/cleaned_pdbs/{name}_clean.pdb"
    resname, atoms = None, []
    for l in open(p):
        if l.startswith("ATOM") and int(l[22:26]) == resseq:
            resname = l[17:20].strip()
            atoms.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return resname, atoms


def nearest_pocket(name, res_atoms):
    best = None
    for pk, sph in C.all_pocket_spheres(name).items():
        md = min(C.dist(s, a) for s in sph for a in res_atoms)
        if best is None or md < best[1]:
            best = (pk, md)
    return best


def raw_for_pocket(name, pk):
    d = C.fpocket_info(name)[pk]
    eo, _, _ = GE.grid_enclosure(name, pk)
    return {"volume": float(d["Volume"]) ** 0.5, "hydrophobicity": float(d["Hydrophobicity score"]),
            "polarity": float(d["Proportion of polar atoms"]), "enclosure": 1 - eo}, float(d["Druggability Score"])


def z_of(m, val):
    mu, sd = Z.stats[m]
    z = (val - mu) / sd
    return -z if m in Z.FLIP else z

def composite_of(raw):
    return sum(Z.WEIGHTS[m] * z_of(m, raw[m]) for m in Z.METRICS) / sum(Z.WEIGHTS.values())


if __name__ == "__main__":
    print(f"{'TF':7} {'residue':9} {'func_pk':8} {'dist(A)':8} {'drugScore':10} "
          f"{'func_comp':10} {'decoy_pk':9} {'decoy_comp':11}")
    for name, (rn, seq) in FUNCTIONAL.items():
        found, atoms = residue_atoms(name, seq)
        if not atoms:
            print(f"{name.split('_')[0]:7} {rn}{seq:<5}  -- residue not resolved in structure, skipped")
            continue
        pk, md = nearest_pocket(name, atoms)
        raw, drug = raw_for_pocket(name, pk)
        fcomp = composite_of(raw)
        decoy_pk = C.PROTEINS[name]["fp_pocket"]
        dcomp = Z.composite(name)
        tag = "" if found == rn else f" (!found {found})"
        print(f"{name.split('_')[0]:7} {rn}{seq:<5}{tag} P{pk:<6} {md:<8.1f} {drug:<10.3f} "
              f"{fcomp:<10.2f} P{decoy_pk:<7} {dcomp:<11.2f}")
