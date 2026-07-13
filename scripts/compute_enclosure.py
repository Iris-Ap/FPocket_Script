"""
compute_enclosure.py  --  Option B: pocket ENCLOSURE computed directly from fpocket
alpha spheres, with NO cross-tool matching (immune to the FOXM1 tool-disagreement issue).

Concept (mirrors DoGSite's lid/hull enclosure): standing inside the pocket, in what
fraction of directions do you hit protein (walled/hull) vs. open space (solvent/lid)?

Algorithm:
  For each alpha-sphere centre c of the selected pocket:
    - take unit directions to all protein atoms within D_NBR of c
    - a look-direction (Fibonacci sphere, N_DIR of them) is "blocked" if some
      protein-atom direction lies within a CONE_DEG half-angle of it
    - burial(c) = blocked / N_DIR
  enclosure = mean burial(c) over the pocket's spheres.  Range 0 (open) .. 1 (buried).
  HIGHER = more enclosed = more druggable  (already oriented; no flip needed).

Validation: prints the DoGSite enclosure of the matched pocket alongside, so the
fpocket-derived values can be checked for consistent ordering on the controls.

Run:  python3 compute_enclosure.py
"""
import os, sys, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C

# ---- parameters (tunable) ----
N_DIR    = 256    # number of look-directions
D_NBR    = 10.0   # Å: neighbourhood to consider protein "walls"
CONE_DEG = 20.0   # deg: half-angle for counting a direction as blocked


def fibonacci_sphere(n):
    pts = []
    phi = math.pi * (3 - math.sqrt(5))
    for i in range(n):
        y = 1 - 2 * (i + 0.5) / n
        r = math.sqrt(max(0.0, 1 - y * y))
        t = phi * i
        pts.append((math.cos(t) * r, y, math.sin(t) * r))
    return pts

DIRS  = fibonacci_sphere(N_DIR)
COS_T = math.cos(math.radians(CONE_DEG))


def burial(center, protein):
    vecs = []
    for p in protein:
        dx, dy, dz = p[0]-center[0], p[1]-center[1], p[2]-center[2]
        d = math.sqrt(dx*dx + dy*dy + dz*dz)
        if 0 < d <= D_NBR:
            vecs.append((dx/d, dy/d, dz/d))
    if not vecs:
        return 0.0
    blocked = 0
    for u in DIRS:
        for v in vecs:
            if u[0]*v[0] + u[1]*v[1] + u[2]*v[2] >= COS_T:
                blocked += 1
                break
    return blocked / N_DIR


def pocket_enclosure(name, pocket):
    spheres = C.alpha_spheres(name, pocket)
    if not spheres:
        return None
    protein = C.protein_atoms(name)
    return sum(burial(c, protein) for c in spheres) / len(spheres)


if __name__ == "__main__":
    print("Enclosure = solid-angle burial from fpocket alpha spheres")
    print(f"params: N_DIR={N_DIR}  D_NBR={D_NBR} A  cone={CONE_DEG} deg\n")
    print(f"{'protein':12} {'type':8} {'fpk':5} {'enclosure':10} "
          f"{'DoGSite_ref':12} {'match(ov)':10}")
    for name in C.ORDER:
        info = C.PROTEINS[name]
        pk = info["fp_pocket"]
        enc = pocket_enclosure(name, pk)
        # DoGSite reference: enclosure of best-overlap pocket (sanity check on controls)
        pid, ov, nsp = C.match_dogsite_pocket(name, pk)
        dref = C.dogsite_desc(name).get(pid, {}).get("enclosure", "-") if pid else "-"
        note = f"{pid}({ov}/{nsp})" if pid else "-"
        print(f"{name:12} {info['type']:8} P{pk:<4} {enc:<10.3f} {dref:12} {note:10}")
