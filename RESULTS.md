# Druggability Results — fpocket vs. FTMap vs. Known Drug Sites

For each protein we compare **where fpocket and FTMap locate a druggable/ligandable site** against the **experimentally known drug-binding site** from the literature. This tests which computational lens correctly identifies real drug sites — especially on transcription-factor DNA-binding domains (DBDs), which are drugged not through classical pockets but at their **DNA-binding surface** (DNA mimicry).

**Scoring conventions**
- fpocket **Druggability Score**: 0–1 (≳0.5 ≈ druggable). *Top pocket = fpocket's #1 ranked; Max = best of all pockets.*
- FTMap: **top consensus-site cluster count** (≳16 ≈ druggable hot spot).
- "Distance to drug site" = distance from the tool's identified site to the known drug's anchor residue / bound-ligand position (as found on literature).

---

## Summary table

| Protein | PDB | Class | Known drug → target site | fpocket | FTMap (clusters) | Correct? |
|---|---|---|---|---|---|---|
| **ERα** | 3ERT | classic druggable ref | tamoxifen → LBD orthosteric pocket | **0.954** ✅ on pocket (1.23 Å) | **21**, on tamoxifen pocket | **both correct** |
| **FOXM1** | 3G73 | TF DBD (liganded) | XST-20 → **His287** (DNA-competitive) | top 0.077 / max 0.181 ❌ | **19**, on DNA surface, **8–13 Å** from His287 | fpocket misses; FTMap right surface, offset |
| **FOXO3** | 2UZK | TF DBD (liganded) | carbenoxolone → **Arg249**; S9 → Arg211/His212/Ser215/Arg222 | **0.856 on a DECOY** (21 Å away); real site **0.001** ❌❌ | **26**, on DNA surface = **S9 site exactly**, **6–10 Å** from Arg249 | fpocket **inverted**; FTMap right surface (matches S9) |
| **FOXA1** | 7VOX | TF DBD (target) | covalent ligands → **C258 (Wing2)** | top 0.027 / max 0.232 ❌ | **30** (strongest of all), on DNA surface (Wing2) | fpocket misses; FTMap on the drugged surface |

**Headline:** fpocket works on the classic target (ERα) but **fails on every TF DBD** — it misses the site (FOXM1, FOXA1) or even scores the *real* drug site as undruggable while hyping a decoy (FOXO3). FTMap consistently lands on the correct **DNA-competitive surface** on all four, though it is offset (~6–13 Å) from the exact anchor residue on the TFs.

---

## ERα (3ERT) — classic druggable reference

- **Known drug:** 4-hydroxytamoxifen (`OHT`), a selective estrogen-receptor modulator, binds the **ligand-binding domain (LBD) orthosteric pocket**. Structure/reference: [RCSB 3ERT](https://www.rcsb.org/structure/3ERT) (Shiau et al., *Cell* 1998).
- **fpocket:** top pocket **Druggability 0.954**; center-to-center distance to the removed drug = **1.23 Å**, with **100 %** of drug atoms within 4 Å → fpocket **blindly rediscovered the real drug pocket.** ✅
- **FTMap:** top consensus site = **21 clusters**; strongest hot-spot residues **Glu353, Arg394, Asp351, His524, Trp383** = the classic tamoxifen-binding pocket. ✅
- **Verdict:** Both methods correct — validates the pipeline on a known druggable target.

## FOXM1 (3G73) — positive control (ligandable forkhead)

- **Known drug:** **XST-20**, identified by structure-based screening, binds the **DBD** at **His287** (plus Arg286, Ser290) — DNA-contacting residues → **DNA-competitive**. Reference: [Zhang et al., *Cell Death Discovery* 2022](https://www.nature.com/articles/s41420-022-01070-w) (PMC9184618).
- **fpocket:** top pocket **0.077**, max **0.181** → **low; misses** the shallow DNA-competitive site.
- **FTMap:** top consensus site = **19 clusters** on the **DNA-binding surface**; strongest residues **Ser306, Trp308, Arg297, Leu259**. His287 (the XST-20 anchor) is **not directly engaged**; Arg286/Ser290 only weakly. Distance from FTMap's peak to His287 = **~8–13 Å**.
- **Verdict:** fpocket misses; FTMap finds the **right surface** but is **offset ~8–13 Å** from the exact XST-20 site.

## FOXO3 (2UZK) — liganded forkhead (the striking case)

- **Known drugs:**
  - **Carbenoxolone (CBX)** binds the DBD at **Arg249** (ionic contact to a DNA-phosphate-gripping residue → DNA mimicry). Reference: [Salcher et al., *Oncogene* 2019](https://www.nature.com/articles/s41388-019-1044-7) (PMC6989399).
  - **S9** (a phenylpyrimidinylguanidine) binds at **Arg211/His212/Ser215/Arg222** (recognition helix). Reference: [Hagenbuchner et al., *eLife* 2019;8:e48876](https://elifesciences.org/articles/48876).
- **fpocket (inverted):** its top-ranked pocket scores **0.856** but sits on a **decoy hydrophobic cleft 21.4 Å from Arg249** (no known ligand); the **actual carbenoxolone pocket** (which lines Arg249 at 4.2 Å) scores **0.001**. → fpocket **rated the real drug site as undruggable and hyped a decoy.** ❌❌
- **FTMap:** top consensus site = **26 clusters** on the DNA surface; strongest residues **Ser215, Arg222, His212, Arg211** — which **match the S9 binding site exactly.** Distance to Arg249 (the CBX site) = **~6–10 Å** (Arg249 itself not engaged).
- **Verdict:** fpocket **inverts** the answer; FTMap lands on the **correct DNA-competitive surface** and precisely recovers the **S9** drug site (offset ~6–10 Å from the CBX site).

## FOXA1 (7VOX) — the target

- **Known drug:** covalent small molecules engage **Cys258 in Wing2** of the forkhead DBD, near the DNA-binding interface and a cancer-mutation-rich region. Reference: [Johnston et al., *Molecular Cell* 2024](https://pubmed.ncbi.nlm.nih.gov/39413792/) (PMC11560529).
- **fpocket:** top pocket **0.027**, max **0.232** (best pocket sits ~12.5 Å from the DNA) → **low; misses** the functional surface.
- **FTMap:** top consensus site = **30 clusters** — the **strongest of all four proteins** — on the **DNA-binding surface**; strongest residues **Ser223, Trp244, His220, Arg219, Ser242, Tyr175, Ile176** (recognition helix + wings). *Note:* the covalent anchor **C258 is beyond the resolved chain** (7VOX chain A = residues 168–252), but FTMap flagged the adjacent **Wing2 residues 242/244** — the same region.
- **Verdict:** fpocket misses; FTMap identifies a **strong hot spot on FOXA1's DNA-binding/Wing2 surface**, matching (at the surface level) the region drugged by covalent ligands in the literature. Notably, FTMap's FOXA1 hot spot overlaps **Ser242**, a recurrent breast-cancer mutation site.

---

## Overall conclusions

1. **Geometry-based scoring (fpocket) fails on transcription-factor DBDs.** It succeeds on the classic target (ERα, 0.954, on the tamoxifen pocket) but on TFs it misses the site (FOXM1, FOXA1) or actively **inverts** it (FOXO3: real drug site 0.001, decoy 0.856).
2. **Interaction-based mapping (FTMap) is far more reliable on TFs.** It lands on the correct **DNA-competitive surface** for all four proteins and, for FOXO3, exactly recovers the S9 drug site.
3. **FTMap is surface-level, not site-exact.** On the TFs its hot spot is offset **~6–13 Å** from the precise known anchor residue (Arg249 in FOXO3, His287 in FOXM1) and does not always engage the exact residue.
4. **FOXA1 is ligandable at its DNA-binding surface.** Its FTMap hot spot is the strongest in the panel (30 clusters), on the Wing2/DNA-recognition surface, coinciding with both the region drugged by covalent ligands and a recurrent cancer mutation (S242).

*Data sources:* raw fpocket output in `raw_FPOCKET_results/`, raw FTMap output in `raw_FTMAP_results/`. Pipeline: `pipeline.sh`.
