"""
drug_enclosure.py  --  Control localization checks: does the selected pocket actually
sit on the bound drug?  (Only the 3 controls have a bound drug.)

For each control it reports, for the fpocket highest-druggability pocket:
  - distance from the pocket's alpha-sphere centroid to the drug centroid
  - fraction of drug atoms within 4 A of the pocket's LINING RESIDUES
    (lining residues = protein residues with an atom within 4 A of an alpha sphere)
  - which fpocket pocket is actually NEAREST the drug (localization sanity check)

Run:  python3 drug_enclosure.py
"""
import os, sys, math, re
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

BASE = C.BASE

def residue_atoms_by_res(name):
    """Protein atoms grouped by (chain,resSeq) from the cleaned PDB."""
    p = f"{BASE}/cleaned_pdbs/{name}_clean.pdb"
    res = {}
    for l in open(p):
        if l.startswith("ATOM"):
            key = (l[21], l[22:26].strip())
            res.setdefault(key, []).append(C._xyz(l))
    return res

def lining_residue_atoms(name, spheres, cutoff=4.0):
    """All atoms of residues that line the given pocket (any atom within cutoff of a sphere)."""
    res = residue_atoms_by_res(name)
    out = []
    for key, atoms in res.items():
        if any(min(C.dist(a, s) for s in spheres) <= cutoff for a in atoms):
            out.extend(atoms)
    return out

if __name__ == "__main__":
    print(f"{'control':7} {'drug':5} {'fpk':5} {'pk->drug(A)':12} "
          f"{'drug_encl(<=4A)':16} {'nearest_pk->drug'}")
    for name in C.ORDER:
        if C.PROTEINS[name]["type"] != "control":
            continue
        drug = C.drug_atoms(name)
        dc = tuple(sum(a[i] for a in drug)/len(drug) for i in range(3))
        pockets = C.all_pocket_spheres(name)
        selpk = C.PROTEINS[name]["fp_pocket"]

        # selected pocket: centroid distance + residue-based drug enclosure
        sp = pockets[selpk]
        cen = tuple(sum(s[i] for s in sp)/len(sp) for i in range(3))
        pk_drug = C.dist(cen, dc)
        lining = lining_residue_atoms(name, sp)
        encl = sum(1 for d in drug if min(C.dist(d, a) for a in lining) <= 4.0) / len(drug) * 100

        # which pocket is nearest the drug (by centroid)
        nearest = min(pockets, key=lambda p:
                      C.dist(tuple(sum(s[i] for s in pockets[p])/len(pockets[p]) for i in range(3)), dc))
        ndist = C.dist(tuple(sum(s[i] for s in pockets[nearest])/len(pockets[nearest]) for i in range(3)), dc)

        print(f"{name.split('_')[0]:7} {C.PROTEINS[name]['drug']:5} P{selpk:<4} "
              f"{pk_drug:<12.2f} {str(round(encl))+'% ('+str(round(encl/100*len(drug)))+'/'+str(len(drug))+')':16} "
              f"P{nearest} ({ndist:.1f} A)")
