# Results

*Placeholders: [X] = table/figure number.*

## 1. fpocket localizes the true drug site more reliably than DoGSiteScorer

To choose a pocket-selection tool, we compared how closely each tool's highest-druggability
pocket sat to the known drug site across the three controls (Table [X]). Both tools were
ranked by their respective druggability scores. fpocket's selected pocket localized to the
bound drug in all three controls — 1.23 Å (ERα), 1.09 Å (HER2), and 6.90 Å (PI3Kα) — with
the pocket enclosing the majority of drug atoms in each case. DoGSiteScorer localized
correctly for only two of three: for ERα and HER2 its top pocket sat on the drug (5.5 and
4.2 Å), but for PI3Kα its highest-drugScore pocket lay 34.1 Å from the drug, enclosing none
of it — a decoy. The true PI3Kα drug pocket ranked below at least eight higher-scoring
DoGSite pockets.

Notably, this failure mode was not unique to DoGSite: fpocket's default detection ranking
also mis-ranked PI3Kα, placing a non-druggable decoy (druggability 0.110) first; the true
drug pocket (druggability 0.899) was only recovered by re-ranking on druggability score.
Both tools also agreed that all three controls are druggable. Given its 3/3 localization
performance — and because fpocket is open-source and scriptable, enabling reproducible
batch analysis where DoGSiteScorer runs as a rate-limited web service — fpocket was adopted
as the pocket-selection tool, with DoGSiteScorer retained for cross-validation.

(The larger PI3Kα distance, 6.90 Å, reflects a centroid offset for a large drug in a large
pocket; atom-level overlap confirmed the pocket fully encloses the drug.)

### Control localization summary

| Control | Drug | fpocket pocket→drug | DoGSite pocket→drug |
|---|---|---|---|
| ERα (3ERT)   | OHT | 1.23 Å | 5.5 Å |
| HER2 (3PP0)  | 03Q | 1.09 Å | 4.2 Å |
| PI3Kα (3HHM) | KWT | 6.90 Å | 34.1 Å (decoy) |

## 2. The geometric rubric ranks all controls above all TFs

Applying the four-metric composite to all six proteins produced a clean separation between
validated drug targets and transcription factors. All three controls scored above all three
TFs. Because the controls define the z-score baseline, they center near zero, while the TFs
fell well below. FOXM1 scored lowest, driven by a small (175 Å³), open, weakly hydrophobic
pocket. At a coarse level, the geometric rubric registers that TF pockets are markedly less
druggable than known drug-target pockets.

| Rank | Protein | Type | Composite |
|---|---|---|---|
| 1 | ERα   | control | +0.82 |
| 2 | HER2  | control | −0.36 |
| 3 | PI3Kα | control | −0.46 |
| 4 | FOXO3 | TF | −0.96 |
| 5 | FOXA1 | TF | −1.72 |
| 6 | FOXM1 | TF | −3.01 |

## 3. Geometry mis-locates druggability: the FOXO3 inversion

FOXO3 was the highest-scoring TF (−0.96) — but this ranking is an artifact of which pocket
the selection rule chose. fpocket's highest-druggability pocket for FOXO3 (druggability
0.856) is a large, hydrophobic, well-enclosed cavity located 21 Å from Arg249, the
functional drug-relevant residue. When the rubric is instead applied to the pocket at the
functional site (fpocket P2, 4.2 Å from Arg249; druggability 0.001), the composite drops
from −0.96 to −3.05 — the lowest of any protein tested, below even FOXM1. The functional
pocket is worse on every axis: smaller (305 vs 674 Å³), far less hydrophobic (7.3 vs 46.4),
more polar (47.6 vs 27.5%), and less enclosed.

This inversion is the central cautionary result: geometric druggability does not merely
mis-estimate the score, it mis-identifies the location, promoting a druggable-looking decoy
over the functional site. Left to select the pocket, geometry reports FOXO3 as having a
druggable pocket (0.856); asked about the site that actually matters, the same geometry
calls it essentially undruggable (0.001).

(FOXM1's functional site, His287, scored −2.62 versus −3.01 for its selected pocket — both
undruggable, with no comparable inversion; FOXA1's functional residue Cys258 is not
resolved in 7VOX. FOXO3 is thus the one clean inversion case.)

| FOXO3 pocket | Location | Volume | Hydrophob | Polarity | drugScore | Composite |
|---|---|---|---|---|---|---|
| Decoy P1 (selected) | 21 Å from Arg249 | 674 | 46.4 | 27.5 | 0.856 | −0.96 |
| Functional P2 | 4.2 Å from Arg249 | 305 | 7.3 | 47.6 | 0.001 | −3.05 |

## 4. The tools do not agree on where TF pockets are

We validated the reimplemented enclosure metric against DoGSiteScorer on pockets both tools
detect in common. Enclosure agreed closely and in identical rank order on the three
controls (this method vs. DoGSite: ERα 0.105/0.06, HER2 0.201/0.07, PI3Kα 0.204/0.19).
However, cross-tool comparison was only possible where the two tools detect the same pocket.
For FOXM1, the fpocket-selected pocket had no overlap with any DoGSite pocket (0 of 17 alpha
spheres within 4 Å of any DoGSite pocket's lining residues); for FOXA1 the overlap was weak
(6/16). On the TFs the two established tools frequently fail to agree even on where the
pocket is — a further indication that these proteins lack well-defined druggable pockets,
and the reason enclosure was computed self-contained rather than transferred from DoGSite.

| Protein | fpocket→DoGSite overlap | Interpretation |
|---|---|---|
| ERα / HER2 / PI3Kα | 44–80 spheres | tools agree strongly |
| FOXO3 | 49/58 spheres | tools agree (both find the decoy) |
| FOXA1 | 6/16 spheres | weak agreement |
| FOXM1 | 0/17 spheres | tools disagree — no shared pocket |
