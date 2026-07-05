# Druggability pipeline

A small pipeline that validates **fpocket** on a protein with a known drug pocket,
as a first step toward assessing the druggability of Forkhead transcription
factors (FOXA1, FOXM1, FOXO3). 

## Validation case
**ERα ligand-binding domain (PDB 3ERT)**, which has the breast-cancer drug
**tamoxifen** (`OHT`) bound in its pocket. We collect its data from PDB (https://www.rcsb.org/structure/3ERT), remove the drug, run fpocket on the bare protein, and check that fpocket independently rediscovers the real drug site. This test was succesful, with the FPocket projection almost exactly lining up with the tamoxifen location. 

## What `pipeline.sh` does

1. **Downloads** `3ERT.pdb` from the RCSB PDB.
2. **Clean** it into `protein_clean.pdb` (bare protein) and `drug_only.pdb`
   (the drug, kept as an "answer key").
3. **Run fpocket** on the bare protein → `protein_clean_out/` (pocket scores).
4. **Check** that fpocket's top pocket sits where the drug was
   (center-to-center distance + fraction of drug atoms inside the pocket).
5. **Prep** `pocket1.pdb` for visualization through Chimera.

Result: Pocket 1 scores **0.954 druggability** and overlaps the drug site
(**1.23 Å** center-to-center, **100%** of drug atoms inside). FPocket is thus validated.

## How to run

```bash
bash pipeline.sh
```

Requires **fpocket** on your PATH (built from https://github.com/Discngine/fpocket).

## Visualize on Chimera (optional)

Open `view_pocket1.cmd` in **UCSF Chimera** (File → Open). It shows the protein
(gray ribbon), tamoxifen (red sticks), and fpocket's Pocket 1 (orange spheres) —
the red drug should sit inside the orange blob.

## Files

Only `pipeline.sh`, `view_pocket1.cmd`, `README.md`, and `.gitignore` are tracked.
All `.pdb` files and `protein_clean_out/` are **generated** by `pipeline.sh`
(see `.gitignore`).
