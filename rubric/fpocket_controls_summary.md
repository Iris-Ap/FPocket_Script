# fpocket results — positive controls (max-druggability pocket)

Pocket selection rule: **the pocket with the highest Druggability Score** in each
structure (not fpocket's default `Score` ranking). Values from
`fpocket_out/<name>/<name>_clean_info.txt`. Inputs = cleaned, single-chain
(chain A), ligand/DNA/water-free structures. Drug coordinates (for the distance check)
from the raw PDBs `fpocket_runs/<PDBID>/<PDBID>.pdb`, chain-A ligand copy only.

## Table 1 — Descriptors of the most-druggable pocket

| Descriptor | ERα (3ERT) | HER2 (3PP0) | PI3Kα (3HHM) |
|---|---|---|---|
| Bound drug (reference) | OHT | 03Q | KWT |
| Selected pocket | P1 | P1 | **P33** |
| **Druggability Score** | **0.954** | **0.526** | **0.899** |
| Volume (Å³) | 857.0 | 1253.8 | 1569.1 |
| Number of α-spheres | 107 | 146 | 156 |
| Total SASA (Å²) | 222.5 | 321.4 | 411.5 |
| Polar SASA (Å²) | 48.4 | 149.6 | 188.7 |
| Apolar SASA (Å²) | 174.1 | 171.8 | 222.8 |
| Mean local hydrophobic density | 57.1 | 39.9 | 29.8 |
| Hydrophobicity score | 58.1 | 28.3 | 14.7 |
| Apolar α-sphere proportion | 0.766 | 0.527 | 0.385 |
| Polarity score | 7 | 14 | 22 |
| Flexibility | 0.176 | 0.193 | 0.295 |

## Table 2 — Localization: does the most-druggable pocket sit on the drug?

| Control | Drug | Selected pocket | Druggability | Distance pocket→drug | Drug atoms enclosed (≤4 Å) |
|---|---|---|---|---|---|
| ERα (3ERT)   | OHT | P1  | 0.954 | 1.23 Å | 69% (20/29) |
| HER2 (3PP0)  | 03Q | P1  | 0.526 | 1.09 Å | 100% (34/34) |
| PI3Kα (3HHM) | KWT | P33 | 0.899 | 6.90 Å | 58% (18/31) |

*Enclosure = fraction of the bound drug's atoms within 4 Å of the pocket's **lining
residues** (residues with an atom within 4 Å of the pocket's alpha spheres) — the same
reference used for the DoGSite table, so the two are directly comparable. All three
pockets enclose the majority of the drug; PI3Kα's larger center-to-center distance is a
centroid offset (large drug/pocket), not a mislocation.*

## Notes

- **Selection rule matters only for PI3Kα.** For ERα and HER2 the most-druggable
  pocket *is* fpocket's Pocket 1. For PI3Kα, fpocket's default Pocket 1 is a
  non-druggable decoy (druggability 0.110, 50.8 Å from the drug); the real drug pocket
  is **P33 (druggability 0.899)**, which fpocket's default `Score` buried at rank 33.
  Sorting by Druggability Score recovers it.
- With this rule, **all three controls select genuinely druggable drug pockets**
  (0.954 / 0.526 / 0.899) that sit on the bound inhibitor (1.1–6.9 Å) — a valid
  druggable baseline.
- The PI3Kα distance (6.90 Å, center-to-center) is looser than ERα/HER2 because KWT is
  large and P33's centroid is offset; atom-overlap confirms P33 is the drug pocket.
