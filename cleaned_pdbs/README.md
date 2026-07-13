# Cleaned PDB inputs — for fpocket & DoGSiteScorer

These are the **exact same files** to feed into **both** fpocket and DoGSiteScorer,
so all pocket/druggability results are directly comparable. Each is a single protein
chain (chain A) with **all DNA, water, drug ligands, ions, and extra protein copies
removed**. Source of each cleaned file: `../fpocket_runs/<PDBID>/protein_clean.pdb`
(copied here 2026-07-12). Raw downloads (with everything still in them) live at
`../fpocket_runs/<PDBID>/<PDBID>.pdb`.

## What was removed from each structure

| File | PDB | Kept | Removed: protein chains | Removed: DNA chains | Removed: ligands / ions | Removed: waters |
|---|---|---|---|---|---|---|
| `ERa_3ERT_clean.pdb`   | 3ERT | chain A | — | — | **OHT** (4-hydroxytamoxifen, 29 atoms) | 79 |
| `HER2_3PP0_clean.pdb`  | 3PP0 | chain A | B | — | **03Q** (inhibitor, 68 atoms) | 190 |
| `PI3Ka_3HHM_clean.pdb` | 3HHM | chain A (p110α) | B (p85α regulatory subunit) | — | **KWT** (inhibitor, 31 atoms) | 115 |
| `FOXA1_7VOX_clean.pdb` | 7VOX | chain A | B, C, H | D, E, F, G | 2 × MG | 446 |
| `FOXM1_3G73_clean.pdb` | 3G73 | chain A | B | C, D | 2 × MG | 163 |
| `FOXO3_2UZK_clean.pdb` | 2UZK | chain A | C | B, D, E, F | — | 209 |

## Final composition (all files)

| File | ATOM records | HETATM | Chains |
|---|---|---|---|
| ERa_3ERT_clean.pdb   | 1962 | 0 | A |
| HER2_3PP0_clean.pdb  | 2290 | 0 | A |
| PI3Ka_3HHM_clean.pdb | 8448 | 0 | A |
| FOXA1_7VOX_clean.pdb | 704  | 0 | A |
| FOXM1_3G73_clean.pdb | 757  | 0 | A |
| FOXO3_2UZK_clean.pdb | 774  | 0 | A |

## Cleaning rule (applied identically to all 6)

1. Keep **chain A only** (drops duplicate protein copies and binding partners).
2. Remove **all nucleic acid** (DNA chains).
3. Remove **all HETATM**: drug ligands, ions (Mg), and water (HOH).
4. Result: a single apo protein chain — no ligand, no DNA, no solvent.

## Important caveat for the drug-site analysis

Because cleaning removes the bound drug, the cleaned files are **apo** (ligand-free).
The **actual drug-site location** for the controls therefore has to come from the
**raw** PDBs, using the bound ligand's coordinates:

| Control | Drug ligand (in raw PDB) |
|---|---|
| ERα (3ERT)   | OHT |
| HER2 (3PP0)  | 03Q |
| PI3Kα (3HHM) | KWT |

Compute the drug's center of mass from `../fpocket_runs/<ID>/<ID>.pdb`, then measure
the distance from each tool's predicted pocket to it. (Raw and cleaned files share the
same coordinate frame — cleaning only deletes atoms, it does not move them.)
