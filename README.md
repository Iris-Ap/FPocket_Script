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
bash pipeline.sh --pdb ID [--drug CODE] [--chain LETTER]
```

Examples:

```bash
bash pipeline.sh --pdb 3ERT --drug OHT   # ERα: has tamoxifen, one chain (validation runs)
bash pipeline.sh --pdb 7VOX --chain A    # FOXA1: no drug, keep chain A (detection only)
```

- `--drug` is optional — omit it if the structure has no drug (the validation step is then skipped).
- `--chain` is optional — use it to keep a single protein chain and drop DNA / extra copies.
- Each run writes to its own folder, `runs/<ID>/`.

Requires **fpocket** on your PATH (built from https://github.com/Discngine/fpocket).

## Visualize on Chimera (optional)

Each run auto-generates a Chimera script at `runs/<ID>/view.cmd`. Open it in
**UCSF Chimera** (File → Open), or from the terminal:

```bash
chimera runs/<ID>/view.cmd
```

It shows the protein (gray ribbon) and fpocket's Pocket 1 (orange spheres); if a
drug was given, the drug appears as red sticks sitting inside the orange pocket.

## Files

Only `pipeline.sh`, `README.md`, and `.gitignore` are tracked. Everything a run
produces (downloaded structure, cleaned files, fpocket output, `view.cmd`) lands
in `runs/<ID>/`, which is gitignored — regenerate it any time with `pipeline.sh`.
All `.pdb` files and `protein_clean_out/` are **generated** by `pipeline.sh`
(see `.gitignore`).
