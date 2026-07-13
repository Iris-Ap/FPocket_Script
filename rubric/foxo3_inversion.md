# FOXO3 — the druggability inversion (decoy vs. functional pocket)

The single cleanest demonstration of the thesis: scoring the *same* protein with the
*same* rubric, but at two different pockets — the geometry-selected decoy vs. the actual
functional drug site — gives opposite conclusions.

## The two pockets

| Pocket | Location | Volume (Å³) | Hydrophobicity | Polarity (%) | Enclosure | fpocket drugScore | **Composite z** |
|---|---|---|---|---|---|---|---|
| **Decoy — P1** (highest-druggability, *selected*) | 21 Å from Arg249 | 674 | 46.4 | 27.5 | 0.64 | 0.856 | **−0.96** |
| **Functional — P2** (nearest Arg249) | 4.2 Å from Arg249 | 305 | 7.3 | 47.6 | 0.50 | 0.001 | **−3.05** |

*(Composites use √volume, following Halgren 2009; see `weights.md`.)*

*Same control baseline and weights as `zscore_composite.py`. Functional pocket = fpocket
pocket nearest Arg249 (the CBX-anchor drug residue). Computed by
`scripts/functional_site_composite.py`.*

## Why the functional site scores far worse — loses on every metric
- **Hydrophobicity 7.3 vs 46.4** — the real site is barely hydrophobic; drugs need greasy pockets.
- **Polarity 47.6 vs 27.5** — the functional site is much more polar (it's a DNA-binding surface built to grip the charged phosphate backbone) → water-friendly, drug-hostile.
- **Volume 305 vs 674** — smaller.
- **Enclosure 0.50 vs 0.64** — more open.

Undruggable on all four axes → composite −3.05.

## The inversion in context

| Protein / pocket | Composite z |
|---|---|
| ERα (control) | +0.82 |
| HER2 (control) | −0.36 |
| PI3Kα (control) | −0.46 |
| **FOXO3 — decoy P1** | **−0.96** |
| FOXA1 | −1.72 |
| FOXM1 | −3.01 |
| **FOXO3 — functional P2** | **−3.05** |

- FOXO3's "best TF" ranking was **entirely the decoy**. Scored at the site that actually
  matters (Arg249), it drops from −0.96 (best TF) to **−3.05 — the worst of all**, below FOXM1.
- The within-protein gap is **~2 SD**, the only difference being whether geometry picks the
  pocket (decoy) or biology does (functional site).

## Interpretation
Geometric druggability does not merely mis-estimate the *score* — it mis-identifies the
*location*, promoting a druggable-looking decoy 21 Å from the real target. Left to choose,
geometry reports FOXO3 as having a druggable pocket (drugScore 0.856); asked about the
functional site, that same geometry calls it essentially undruggable (0.001). This is why
geometric pocket-druggability cannot be trusted as a standalone readout for transcription
factors.

## Note on the other TFs
- **FOXM1 (His287):** functional residue *is* resolved in 3G73, so an analogous
  functional-site comparison is possible if needed.
- **FOXA1 (Cys258):** the functional residue is **not resolved** in structure 7VOX
  (chain A spans residues 168–252; Cys258 in Wing2 is absent), so FOXA1 cannot receive the
  same functional-site treatment from this structure.
