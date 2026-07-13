"""
grid_enclosure.py -- Pocket enclosure via the DoGSite / Volkamer (2012) grid definition,
reimplemented on the fpocket-selected pocket so EVERY protein gets a value with NO
cross-tool matching.

enclosure = ratio of solvent-exposed "lid" grid points to pocket-boundary "hull" points
(Volkamer et al., Bioinformatics 2012, 28(15):2074).

Reimplementation (v2 -- with solvent probe + bulk-connectivity, to match DoGSite's
solvent-accessible definition):
  1. Voxelize a grid (GRID A) over the fpocket pocket's alpha spheres.
  2. POCKET voxel  = inside the alpha-sphere union AND not inside protein (cavity void).
  3. PROBE-ACCESSIBLE voxel = not pocket AND not within (VDW+PROBE) of protein
     (a 1.4 A solvent probe cannot enter gaps narrower than the probe -> narrow
     openings are closed off, exactly as in a solvent-accessible surface).
  4. BULK = the probe-accessible component(s) that reach the box boundary (true outside
     solvent); internal voids do NOT count.
  5. HULL = pocket voxels with >=1 non-pocket neighbour (boundary shell).
     LID  = hull voxels with >=1 BULK neighbour (boundary open to outside solvent).
  6. enclosure_open = LID / HULL       (DoGSite direction: higher = more OPEN)
     burial         = 1 - enclosure_open (rubric: higher = more ENCLOSED = druggable)

Run: python3 grid_enclosure.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import common as C
import numpy as np
try:
    from scipy import ndimage
    HAVE_SCIPY = True
except ImportError:
    HAVE_SCIPY = False

GRID   = 1.0   # A grid spacing
VDW    = 1.8   # A united-atom protein radius
PROBE  = 1.4   # A solvent probe radius
MARGIN = 5.0   # A padding beyond the pocket


def alpha_spheres_radii(name, pocket):
    fn = f"{C.BASE}/fpocket_out/{name}/pockets/pocket{pocket}_vert.pqr"
    out = []
    for l in open(fn):
        if l.startswith(("ATOM", "HETATM")):
            f = l.split()
            out.append((float(f[-5]), float(f[-4]), float(f[-3]), float(f[-1])))
    return np.array(out)


def _within(P, atoms, r, chunk=300):
    """Boolean mask: grid point within r of any atom (chunked over atoms)."""
    m = np.zeros(P.shape[0], bool)
    r2 = r * r
    for i in range(0, len(atoms), chunk):
        a = atoms[i:i+chunk]
        m |= (((P[:, None, :] - a[None, :, :]) ** 2).sum(-1) <= r2).any(1)
    return m


def grid_enclosure(name, pocket):
    S = alpha_spheres_radii(name, pocket)
    prot = np.array(C.protein_atoms(name))
    cen, rad = S[:, :3], S[:, 3]

    lo = cen.min(0) - rad.max() - MARGIN
    hi = cen.max(0) + rad.max() + MARGIN
    xs = np.arange(lo[0], hi[0], GRID); ys = np.arange(lo[1], hi[1], GRID); zs = np.arange(lo[2], hi[2], GRID)
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    shape = gx.shape
    P = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], 1)

    prot = prot[((prot >= lo - (VDW+PROBE)) & (prot <= hi + (VDW+PROBE))).all(1)]

    in_sphere = (((P[:, None, :] - cen[None, :, :]) ** 2).sum(-1) <= (rad[None, :] ** 2)).any(1)
    in_prot   = _within(P, prot, VDW)              # protein interior
    in_excl   = _within(P, prot, VDW + PROBE)      # probe-inaccessible shell

    POCK = (in_sphere & ~in_prot).reshape(shape)
    PROBE_ACC = (~POCK.ravel() & ~in_excl).reshape(shape)   # probe-accessible solvent

    # BULK = probe-accessible components that touch the box boundary
    if HAVE_SCIPY:
        lbl, _ = ndimage.label(PROBE_ACC, structure=ndimage.generate_binary_structure(3, 1))
        face = np.unique(np.concatenate([
            lbl[0].ravel(), lbl[-1].ravel(), lbl[:, 0].ravel(),
            lbl[:, -1].ravel(), lbl[:, :, 0].ravel(), lbl[:, :, -1].ravel()]))
        face = face[face != 0]
        BULK = np.isin(lbl, face)
    else:
        BULK = PROBE_ACC   # fallback: no connectivity (probe only)

    Pk = np.pad(POCK, 1, constant_values=False)
    Bk = np.pad(BULK, 1, constant_values=True)     # outside box = bulk solvent
    core = Pk[1:-1, 1:-1, 1:-1]
    not_pock_nb = np.zeros_like(core)
    bulk_nb = np.zeros_like(core)
    for dx, dy, dz in [(1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)]:
        sl = (slice(1+dx, Pk.shape[0]-1+dx), slice(1+dy, Pk.shape[1]-1+dy), slice(1+dz, Pk.shape[2]-1+dz))
        not_pock_nb |= ~Pk[sl]
        bulk_nb |= Bk[sl]

    HULL = core & not_pock_nb
    LID = HULL & bulk_nb
    hull, lid = int(HULL.sum()), int(LID.sum())
    return (lid / hull if hull else float("nan")), hull, lid


if __name__ == "__main__":
    print("Grid enclosure v2 (Volkamer 2012 def + 1.4 A solvent probe + bulk connectivity)")
    print(f"params: grid={GRID} vdw={VDW} probe={PROBE} margin={MARGIN} | scipy={HAVE_SCIPY}\n")
    print(f"{'protein':12} {'fpk':5} {'enc_open':9} {'burial':8} {'DoGSite_ref':12} {'match'}")
    for name in C.ORDER:
        pk = C.PROTEINS[name]["fp_pocket"]
        eo, hull, lid = grid_enclosure(name, pk)
        pid, ov, nsp = C.match_dogsite_pocket(name, pk)
        dref = C.dogsite_desc(name).get(pid, {}).get("enclosure", "-") if pid else "-"
        note = f"{pid}({ov}/{nsp})" if pid else "-"
        print(f"{name:12} P{pk:<4} {eo:<9.3f} {1-eo:<8.3f} {dref:12} {note}")
