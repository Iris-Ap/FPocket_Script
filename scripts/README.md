# Analysis scripts — TF druggability rubric

Reproducible code for every number in the rubric. All scripts import `common.py`
(shared config + PDB/fpocket/DoGSite loaders) and run against the folders in `~/SRA_iris`:
`cleaned_pdbs/`, `fpocket_out/`, `dogsite_out/`, `fpocket_runs/` (raw PDBs with drug).

Run any script from this folder, e.g. `python3 compute_enclosure.py`.

| Script | What it does |
|---|---|
| `common.py` | Shared protein registry (which fpocket pocket is selected per protein, bound-drug codes) + loaders for alpha spheres, protein atoms, fpocket descriptors, bound-drug atoms, DoGSite descriptors/residues, and the fpocket↔DoGSite pocket matcher. |
| `pull_rubric_data.py` | Extracts the rubric metrics (Volume, Hydrophobicity, Polarity) from each protein's highest-druggability fpocket pocket. |
| `compute_enclosure.py` | **Option B** — computes pocket enclosure (solid-angle burial) directly from fpocket alpha spheres, no cross-tool matching. Prints DoGSite enclosure alongside for validation. |
| `match_pockets.py` | Matches each fpocket-selected pocket to a DoGSite pocket (alpha-sphere/residue overlap) and reads DoGSite enclosure/depth/drugScore. Overlap = 0 means the tools disagree on the pocket (e.g. FOXM1). |
| `drug_enclosure.py` | Control localization: does the selected pocket sit on the bound drug? Distance + % drug atoms enclosed + which pocket is actually nearest the drug. |

## Pocket selection rule
Rubric pocket = the fpocket pocket with the **highest Druggability Score** (not fpocket's
default `Score` rank). Recorded per protein in `common.PROTEINS[...]["fp_pocket"]`.

## Metric orientation (for z-scoring)
- Volume ↑ (higher = more druggable)
- Hydrophobicity ↑
- Polarity ↓ (flip — more polar = less druggable)
- Enclosure ↑ (as defined in `compute_enclosure.py`: higher burial = more enclosed = more druggable)
