# Methodology

*Placeholders: [CITE] = reference needed; [X] = table/figure number. Suggested refs:
Le Guilloux 2009 (fpocket); Volkamer 2012 (DoGSiteScorer); Halgren 2009 (Dscore);
Cheng 2007 (maximal affinity); Schmidtke & Barril 2010 (fpocket druggability);
Lee & Richards 1971 (solvent probe).*

## 1. Study design and tool selection

Transcription factors (TFs) are widely regarded as difficult to drug, in large part
because they often lack the deep, well-defined binding pockets that characterize classical
drug targets. To evaluate whether geometry-based pocket-detection methods can meaningfully
assess the druggability of forkhead-box (FOX) transcription factors, we built a scoring
rubric that produces a relative composite druggability ranking, calibrated against a set of
validated drug targets used as positive controls.

The rubric combines four pocket features: hydrophobicity, polarity, pocket volume, and
enclosure. These descriptors are established determinants of druggability — hydrophobic,
enclosed pockets favor high-affinity ligand binding by promoting desolvation and shielding
the ligand from water, whereas polar, solvent-exposed surfaces resist small-molecule
binding [CITE: Cheng 2007; Halgren 2009; Schmidtke & Barril 2010].

**Pocket detection.** To score a protein, we must first identify and select a pocket on it.
We use two established pocket-detection tools: fpocket, which identifies cavities using
Voronoi tessellation and alpha spheres [CITE: Le Guilloux 2009], and DoGSiteScorer, which
detects pockets on a grid using a difference-of-Gaussians filter [CITE: Volkamer 2012].
Both take a cleaned protein structure, return a ranked list of pockets, and report
per-pocket descriptors (volume, hydrophobicity, surface area, enclosure, and a druggability
score, among others). We selected each protein's highest-druggability pocket for analysis.
Notably, fpocket's default ranking is by an internal detection score rather than
druggability, so we re-ranked pockets by the reported druggability score.

**Choosing which tool to trust.** Because fpocket and DoGSiteScorer may return different
top-ranked pockets for the same protein, we first determined which tool more reliably
identifies the true druggable pocket. We used three positive controls — ERα (3ERT), HER2
(3PP0), and PI3Kα (3HHM) — all validated drug targets whose ligand-binding sites are known
from co-crystallized inhibitors. Each structure contains its bound drug (3ERT→OHT,
4-hydroxytamoxifen; 3PP0→03Q; 3HHM→KWT); we used the center of mass of each bound drug as
ground truth for the true druggable-site location. We then removed the drug, ran the
ligand-free structure through both tools, and measured the distance (Å) between each tool's
selected pocket and the drug's center of mass. The tool whose selected pocket most
consistently localizes to the true drug site across the three controls was adopted as our
pocket-selection method for the TFs.

## 2. Structure preparation

All structures were obtained from the Protein Data Bank and prepared identically. For each
protein we retained chain A only and removed all non-protein components: nucleic acids
(DNA), water, ions, and any bound ligand. Chain A was kept because it is the first-listed
and typically most complete protein copy in the asymmetric unit; restricting to a single
chain removes redundant crystallographic copies and binding partners so that each protein's
intrinsic pockets are assessed in isolation. For the FOX transcription factors, the
deposited structures are protein–DNA complexes; DNA was removed so that pocket detection
reflects the protein surface rather than the protein–DNA interface, and so that the TFs and
controls are treated by an identical cleaning procedure. The bound drug was removed from
each control prior to pocket detection so that the tools predict pockets de novo, without
being biased by the ligand; the removed drug's coordinates were retained separately to
define the ground-truth drug-site location (Section 1). The exact components removed from
each control are listed in Table [X].

## 3. Rubric metrics

The rubric scores a single pocket using four descriptors. Three are taken directly from
fpocket's output for the selected pocket:

- **Pocket volume** (Å³), entered as √volume following Halgren's druggability model, which
  uses √(pocket size); the square-root reflects the diminishing return of increasing pocket
  size beyond the minimum required to accommodate a drug [CITE: Halgren 2009].
- **Hydrophobicity** — fpocket's residue-based hydrophobicity score, reflecting the
  greasiness of the amino acids lining the pocket.
- **Polarity** — the proportion of polar atoms lining the pocket. This is a distinct,
  atom-level property from hydrophobicity (a residue-level property): a pocket lined by
  hydrophobic residues can still expose polar backbone or side-chain atoms, so the two
  descriptors capture related but non-identical chemistry.

The fourth descriptor, enclosure, is not reported by fpocket. Rather than take it from
DoGSiteScorer — which would require matching fpocket's pocket to a corresponding DoGSite
pocket, and is not always possible (Results) — we computed enclosure directly on the
fpocket-selected pocket, reimplementing the grid-based definition of Volkamer et al.
[CITE: Volkamer 2012]. The pocket cavity was voxelized on a 1 Å grid from fpocket's alpha
spheres. Each voxel was classified as pocket (inside the alpha-sphere cavity, not
overlapping protein), protein (within 1.8 Å of a protein atom), or solvent. Solvent
accessibility was determined with a 1.4 Å probe [CITE: Lee & Richards 1971], excluding
voxels within 3.2 Å (van der Waals + probe) of protein, and retaining only the solvent
region connected to the exterior. Enclosure was computed as the ratio of solvent-exposed
("lid") boundary voxels to total pocket-boundary ("hull") voxels; we report
burial = 1 − (lid/hull), so that higher values indicate more enclosed pockets. This
reimplementation reproduced DoGSiteScorer's enclosure on pockets both tools detect in
common (Results).

## 4. Composite scoring

Each descriptor was standardized as a z-score relative to the three control proteins, which
define the baseline distribution of a validated druggable pocket:

    z_i = (x_i − μ_control) / σ_control

where μ and σ are the mean and sample standard deviation of that descriptor across the
three controls. Descriptors were oriented so that higher z corresponds to greater
druggability: volume, hydrophobicity, and enclosure (burial) are used directly, while
polarity is sign-flipped, since higher polarity indicates lower druggability.

The four z-scores were combined into a single composite as a weighted average, with weights
reflecting the relative importance of each descriptor in established druggability models —
enclosure 90, hydrophobicity 90, polarity 60, volume 40. This ordering mirrors the
coefficient hierarchy of Halgren's Dscore (Dscore = 0.094·√n + 0.60·enclosure −
0.324·hydrophilicity), in which enclosure carries the largest weight, the
polarity/hydrophilicity penalty roughly half of that, and pocket size the smallest
[CITE: Halgren 2009; Cheng 2007]. A higher composite indicates a pocket more geometrically
consistent with known druggable sites; because the controls define the baseline, they
center near zero, and TF pockets are interpreted relative to this control distribution.
