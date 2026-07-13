"""
Shared config + loaders for the SRA_iris druggability-rubric pipeline.

Folder layout assumed:
  ~/SRA_iris/cleaned_pdbs/<name>_clean.pdb        cleaned chain-A inputs
  ~/SRA_iris/fpocket_out/<name>/                  fpocket output (run on cleaned files)
  ~/SRA_iris/dogsite_out/<name>/                  DoGSiteScorer output
  ~/SRA_iris/fpocket_runs/<PDBID>/<PDBID>.pdb     raw downloads (with bound drug/DNA)

<name> = "<Protein>_<PDBID>" e.g. ERa_3ERT, FOXA1_7VOX.
"""
import os, math, glob, re

BASE = os.path.expanduser("~/SRA_iris")

# Protein registry. fp_pocket = fpocket pocket with the HIGHEST Druggability Score
# (the rubric's selection rule). drug = bound-ligand HETATM code in the raw PDB.
PROTEINS = {
    "FOXA1_7VOX": {"pdbid": "7VOX", "type": "TF",      "fp_pocket": 6,  "drug": None},
    "FOXM1_3G73": {"pdbid": "3G73", "type": "TF",      "fp_pocket": 2,  "drug": None},
    "FOXO3_2UZK": {"pdbid": "2UZK", "type": "TF",      "fp_pocket": 1,  "drug": None},
    "ERa_3ERT":   {"pdbid": "3ERT", "type": "control", "fp_pocket": 1,  "drug": "OHT"},
    "HER2_3PP0":  {"pdbid": "3PP0", "type": "control", "fp_pocket": 1,  "drug": "03Q"},
    "PI3Ka_3HHM": {"pdbid": "3HHM", "type": "control", "fp_pocket": 33, "drug": "KWT"},
}
ORDER = ["FOXA1_7VOX", "FOXM1_3G73", "FOXO3_2UZK", "ERa_3ERT", "HER2_3PP0", "PI3Ka_3HHM"]


def _xyz(l):
    return (float(l[30:38]), float(l[38:46]), float(l[46:54]))

def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


# ---- fpocket ----
def fpocket_out_pdb(name):
    return f"{BASE}/fpocket_out/{name}/{name}_clean_out.pdb"

def alpha_spheres(name, pocket):
    """Alpha-sphere centres of one fpocket pocket (STP atoms with resSeq == pocket)."""
    return [_xyz(l) for l in open(fpocket_out_pdb(name))
            if l.startswith("HETATM") and l[17:20].strip() == "STP" and int(l[22:26]) == pocket]

def all_pocket_spheres(name):
    """{pocket# : [sphere xyz, ...]} for every fpocket pocket."""
    d = {}
    for l in open(fpocket_out_pdb(name)):
        if l.startswith("HETATM") and l[17:20].strip() == "STP":
            d.setdefault(int(l[22:26]), []).append(_xyz(l))
    return d

def protein_atoms(name):
    """Protein atom coords from the cleaned PDB."""
    p = f"{BASE}/cleaned_pdbs/{name}_clean.pdb"
    return [_xyz(l) for l in open(p) if l.startswith("ATOM")]

def fpocket_info(name):
    """Parse fpocket info.txt -> {pocket# : {descriptor: value}}."""
    info = f"{BASE}/fpocket_out/{name}/{name}_clean_info.txt"
    b, cur = {}, None
    for line in open(info):
        m = re.match(r'^Pocket (\d+) :', line)
        if m:
            cur = int(m.group(1)); b[cur] = {}; continue
        if cur and ':' in line:
            k, _, v = line.partition(':'); b[cur][k.strip()] = v.strip()
    return b

def highest_druggability_pocket(name):
    """Pocket number with the max Druggability Score (rubric selection rule)."""
    b = fpocket_info(name)
    return max(b, key=lambda p: float(b[p]["Druggability Score"]))


# ---- bound drug (controls) ----
def drug_atoms(name):
    """Bound-drug heavy atoms (chain A copy) from the raw PDB; None if apo/TF."""
    lig = PROTEINS[name]["drug"]
    if not lig:
        return None
    pid = PROTEINS[name]["pdbid"]
    raw = f"{BASE}/fpocket_runs/{pid}/{pid}.pdb"
    return [_xyz(l) for l in open(raw)
            if l.startswith("HETATM") and l[17:20].strip() == lig and l[21] == "A"]


# ---- DoGSite ----
def dogsite_desc(name):
    """Parse DoGSite _desc.txt -> {P_x : {column: value}}."""
    desc = glob.glob(f"{BASE}/dogsite_out/{name}/*_desc.txt")[0]
    rows, hdr = {}, None
    for l in open(desc):
        p = l.rstrip("\n").split("\t")
        if hdr is None:
            hdr = p; continue
        rows[p[0]] = {hdr[i]: p[i] for i in range(len(hdr))}
    return rows

def dogsite_pocket_residues(name):
    """DoGSite per-pocket lining-residue atom coords -> {P_x : [xyz, ...]}."""
    d = {}
    for fn in glob.glob(f"{BASE}/dogsite_out/{name}/residues/*_res.pdb"):
        pid = re.search(r'_(P_\d+)_res', fn).group(1)
        d[pid] = [_xyz(l) for l in open(fn) if l.startswith(("ATOM", "HETATM"))]
    return d

def match_dogsite_pocket(name, pocket, cutoff=4.0):
    """
    Best-overlapping DoGSite pocket for an fpocket pocket, by alpha-sphere containment.
    Returns (P_x, overlap_count, n_spheres). overlap_count == 0 => tools disagree
    (no shared pocket), e.g. FOXM1.
    """
    spheres = alpha_spheres(name, pocket)
    dogs = dogsite_pocket_residues(name)
    best = (None, -1)
    for pid, res in dogs.items():
        ov = sum(1 for s in spheres if min(dist(s, a) for a in res) <= cutoff)
        if ov > best[1]:
            best = (pid, ov)
    return best[0], best[1], len(spheres)
