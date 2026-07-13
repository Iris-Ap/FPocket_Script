# Rubric weights & justification

Composite = weighted sum of oriented z-scores. Weights reflect the published druggability
hierarchy (enclosure ≈ hydrophobicity > polarity > volume).

| Metric | Weight | Orientation | Justification | Key refs |
|---|---|---|---|---|
| **Enclosure** | 90 | higher = druggable | Burial concentrates hydrophobic contacts and shields ligand from water — largest positive coefficient in Halgren's SiteMap druggability score. | Halgren 2009; Cheng 2007 |
| **Hydrophobicity** | 90 | higher = druggable | Hydrophobic/desolvation effect sets maximal achievable affinity; top feature in fpocket's druggability model. | Cheng 2007; Schmidtke & Barril 2010 |
| **Polarity** | 60 | flip (higher = less druggable) | Polar/charged surfaces bind water and resist hydrophobic drugs; hydrophilicity penalty is ~half the enclosure weight in Halgren's model. | Halgren 2009; Cheng 2007 |
| **Pocket volume** | 40 | higher = druggable | Pocket must be large enough for a drug, but above threshold size is a weak discriminator (smallest coefficient in Halgren's model). Entered as **√volume**, matching Halgren's √n term (diminishing returns of size). | Halgren 2009 |

**Anchor:** Halgren (2009) Dscore = 0.094·√size + 0.60·enclosure − 0.324·hydrophilicity —
the same hierarchy (enclosure > polarity > size) our 90/60/40 weights reproduce.

**References**
- Cheng et al. *Nat. Biotechnol.* 2007, 25, 71.
- Halgren. *J. Chem. Inf. Model.* 2009, 49, 377.
- Schmidtke & Barril. *J. Med. Chem.* 2010, 53, 5858.

**Caveat:** the relative *ordering* of the weights is literature-grounded; the exact
values (90/90/60/40) are a reasonable operationalization of that hierarchy, not lifted
verbatim. Ranking should be checked for stability under reweighting (sensitivity test).
