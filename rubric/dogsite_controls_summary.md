# DoGSiteScorer results — positive controls (highest-drugScore pocket)

Pocket selection rule: **the pocket with the highest `drugScore`** in each structure.
Values from the DoGSite `_desc.txt` files in `dogsite_out/<protein>/`. Inputs = the same
cleaned, chain-A, ligand/DNA/water-free files used for fpocket. Run settings: chain A,
"Properties and druggability", de novo (no ligand). Drug coordinates for the distance
check from `fpocket_runs/<PDBID>/<PDBID>.pdb`, chain-A ligand copy only.

## Table 1 — Descriptors of the highest-drugScore pocket

| Descriptor | ERα (P_0) | HER2 (P_0) | PI3Kα (P_4) |
|---|---|---|---|
| **drugScore** | **0.822** | **0.808** | **0.872** |
| simpleScore | 0.71 | 0.63 | 0.44 |
| Volume (Å³) | 1217.9 | 1260.7 | 746.8 |
| Surface (Å²) | 1460.9 | 1404.5 | 972.4 |
| Depth (Å) | 20.8 | 32.5 | 28.2 |
| Enclosure | 0.06 | 0.07 | 0.09 |
| Surface/Volume | 1.20 | 1.11 | 1.30 |
| Ellipsoid c/a | 0.10 | 0.08 | 0.05 |
| Ellipsoid b/a | 0.25 | 0.28 | 0.20 |
| Site atoms (lining) | 256 | 235 | 195 |
| H-bond acceptors | 59 | 67 | 53 |
| H-bond donors | 15 | 26 | 21 |
| Hydrophobic interactions | 92 | 51 | 26 |
| Hydrophobicity | 0.55 | 0.35 | 0.26 |
| Metal | 0 | 0 | 0 |
| C atoms | 192 | 165 | 136 |
| N atoms | 30 | 30 | 24 |
| O atoms | 30 | 38 | 34 |
| S atoms | 4 | 2 | 1 |
| Other atoms | 0 | 0 | 0 |
| Negative AA (frac) | 0.12 | 0.08 | 0.17 |
| Positive AA (frac) | 0.12 | 0.17 | 0.10 |
| Polar AA (frac) | 0.10 | 0.29 | 0.28 |
| Apolar AA (frac) | 0.67 | 0.46 | 0.45 |

## Table 2 — Distance from the highest-drugScore pocket to the drug

| Protein | Pocket | drugScore | Distance pocket→drug | Drug atoms enclosed (≤4 Å) |
|---|---|---|---|---|
| ERα (3ERT)   | P_0 | 0.822 | 5.5 Å  | 69% (20/29) |
| HER2 (3PP0)  | P_0 | 0.808 | 4.2 Å  | 100% (34/34) |
| PI3Kα (3HHM) | P_4 | 0.872 | **34.1 Å** | **0% (0/31)** |

*Distance = pocket-lining-residue centroid → drug centroid. Enclosure = fraction of drug
atoms within 4 Å of the pocket's lining residues (same reference as the fpocket table).
For ERα/HER2 the highest-drugScore pocket sits on the drug; for **PI3Kα the highest-drugScore
pocket (P_4) is a decoy 34 Å away enclosing none of the drug** — the actual KWT pocket is
P_9 (drugScore 0.703).*
