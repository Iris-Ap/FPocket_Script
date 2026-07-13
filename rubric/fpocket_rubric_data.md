# Druggability rubric — raw data

**Rubric metrics: pocket Volume, Hydrophobicity, Polarity, Enclosure.** (Druggability
Score is used *solely* to select which pocket to analyze — the highest-druggability
pocket in each structure — and is NOT a rubric criterion.)

Pocket selected per protein via fpocket (highest Druggability Score); all geometry runs
on the cleaned, chain-A, DNA/water/drug-free files (`cleaned_pdbs/<name>_clean.pdb`).
Volume/Hydrophobicity/Polarity are raw fpocket descriptors (`fpocket_out/<name>/`).
**Enclosure** is computed on that same fpocket pocket by `scripts/grid_enclosure.py`
(grid hull/lid method per Volkamer 2012 + 1.4 Å solvent probe; see
`rubric/enclosure_methods.md`). Values pre-standardization.

## Raw values

| Protein | Type | Pocket | Volume (Å³) | Hydrophobicity | Polarity (% polar atoms) | Enclosure (burial) |
|---|---|---|---|---|---|---|
| FOXA1 | TF | P6 | 172.6 | 44.0 | 21.4 | 0.543 |
| FOXM1 | TF | P2 | 175.1 | 12.6 | 35.3 | 0.461 |
| FOXO3 | TF | P1 | 673.9 | 46.4 | 27.5 | 0.636 |
| ERα | control | P1 | 857.0 | 58.1 | 25.8 | 0.895 |
| HER2 | control | P1 | 1253.8 | 28.3 | 39.7 | 0.799 |
| PI3Kα | control | P33 | 1569.1 | 14.7 | 39.6 | 0.796 |

## Oriented z-scores (√volume) and composite

z = (value − control mean) / control SD, oriented so higher = more druggable (polarity
flipped). Volume enters as **√volume** (following Halgren 2009). Composite = weighted
average of z-scores (weights: enclosure 90, hydrophobicity 90, polarity 60, volume 40).
Baseline = the 3 controls. From `scripts/zscore_composite.py`.

| Protein | Type | Volume | Hydrophob | Polarity | Enclosure | **Composite** |
|---|---|---|---|---|---|---|
| FOXA1 | TF | −4.16 | 0.46 | 1.70 | −5.10 | **−1.72** |
| FOXM1 | TF | −4.14 | −0.95 | −0.03 | −6.55 | **−3.01** |
| FOXO3 | TF | −1.69 | 0.57 | 0.94 | −3.43 | **−0.96** |
| ERα | control | −1.06 | 1.10 | 1.15 | 1.15 | **+0.82** |
| HER2 | control | 0.12 | −0.24 | −0.59 | −0.55 | **−0.36** |
| PI3Kα | control | 0.93 | −0.86 | −0.57 | −0.61 | **−0.46** |

Ranked (most → least druggable): ERα +0.82 > HER2 −0.36 > PI3Kα −0.46 > FOXO3 −0.96 >
FOXA1 −1.72 > FOXM1 −3.01. **All three controls rank above all three TFs.**

## Metric mapping & orientation (for z-scoring)

| Rubric criterion | Source field / method | Orientation (higher = more druggable?) |
|---|---|---|
| **Volume** | fpocket `Volume`, entered as **√Volume** (Halgren 2009) | ↑ yes |
| **Hydrophobicity** | fpocket `Hydrophobicity score` | ↑ yes |
| **Polarity** | fpocket `Proportion of polar atoms` | ↓ **flip** (more polar = less druggable) |
| **Enclosure** | `grid_enclosure.py` burial = 1 − lid/hull | ↑ yes (higher burial = more enclosed) |

*Volume enters as √Volume, matching Halgren's Dscore (0.094·√n + 0.60·e − 0.324·p); the
√ compresses the top end (diminishing returns of size). Ranking is unchanged vs. linear
volume.*
*Polarity uses `Proportion of polar atoms` (%) not fpocket's `Polarity score` (a raw
count), because the raw count scales with pocket size and would double-count with Volume.*
*Enclosure "burial" = 1 − (lid ÷ hull grid points); already oriented so higher = more
druggable (no flip).*

## Notes

- **Enclosure method & validation:** grid-based hull/lid enclosure reimplemented per
  Volkamer et al. (2012) with a 1.4 Å solvent probe, computed on the fpocket pocket so
  every protein gets a value with no cross-tool matching. Validated against DoGSiteScorer
  on the pockets both tools detect in common (ERα, HER2, PI3Kα — close values + identical
  rank order). FOXO3 is a known cross-tool exception (mine 0.36 open vs DoGSite 0.08);
  FOXA1/FOXM1 have no shared DoGSite pocket to compare. See `rubric/enclosure_methods.md`.
- **FOXO3 caveat:** its selected (highest-druggability) pocket P1 is the **decoy** ~21 Å
  from the real drug site (Arg249) — an artifact of the highest-druggability selection
  rule. Flag it when interpreting FOXO3's rubric score.
