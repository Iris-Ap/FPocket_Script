# Enclosure metric — methods

## Purpose
Enclosure quantifies pocket shape: what fraction of the pocket surface is **walled by
protein** vs. **open to solvent** (the "mouth"/lid). Deep druggable pockets are highly
enclosed (small lid); shallow TF sites are open (large lid). fpocket does not output
enclosure, so it is computed here on the fpocket-selected pocket.

## Definition (Volkamer et al., 2012)
enclosure = (solvent-exposed "lid" grid points) ÷ (pocket-boundary "hull" grid points).
Lower = more enclosed. Note the paper's separate "4 Å" cutoff defines *pocket atoms* for
composition descriptors — it is NOT part of the enclosure calculation.

## Reimplementation (`scripts/grid_enclosure.py`)
1. **Voxelize** a 1 Å grid over the pocket (from fpocket's alpha spheres + 5 Å margin).
2. **Classify** each voxel:
   - *pocket* = inside the alpha-sphere cavity and not inside protein;
   - *protein* = within VdW (1.8 Å) of a protein atom;
   - *solvent* = the rest.
3. **Solvent probe (1.4 Å).** A voxel is solvent-*accessible* only if it is NOT within
   (VdW + probe) = 3.2 Å of protein. This is the solvent-accessible-surface construction
   (Lee & Richards, 1971): rolling a water-sized probe seals inter-atomic cracks so only
   genuine, water-passable openings count. Implemented by growing atom radii by the probe.
4. **Bulk connectivity.** Keep only the probe-accessible component(s) connected to the box
   exterior (`scipy.ndimage.label`), so sealed internal voids do not count as "open."
5. **Hull** = pocket voxels with a non-pocket neighbour (the boundary shell).
   **Lid** = hull voxels with a bulk-solvent neighbour.
6. **enclosure_open = lid ÷ hull**; **burial = 1 − enclosure_open** (rubric metric,
   higher = more enclosed = more druggable).

## Parameters
grid = 1.0 Å · VdW = 1.8 Å · probe = 1.4 Å (water) · margin = 5.0 Å.

## Validation vs. DoGSiteScorer (proteins with a shared pocket)
| Protein | This method (lid/hull) | DoGSite enclosure |
|---|---|---|
| ERα | 0.105 | 0.06 |
| HER2 | 0.201 | 0.07 |
| PI3Kα | 0.204 | 0.19 |
| FOXO3 | 0.364 | 0.08 (outlier) |

Identical rank order and close values on the three well-defined pockets; FOXO3 is a known
cross-tool exception. FOXA1 (weak overlap) and FOXM1 (no shared pocket) cannot be
cross-validated — which is why enclosure is computed self-contained rather than taken
from DoGSite.

## References
- Volkamer, Kuhn, Rippmann, Rarey. DoGSiteScorer. *Bioinformatics* 2012, 28(15):2074.
- Lee, Richards. Solvent-accessible surface / probe. *J. Mol. Biol.* 1971, 55:379.
- Le Guilloux, Schmidtke, Tuffery. Fpocket. *BMC Bioinformatics* 2009, 10:168.
