#!/usr/bin/env bash
#
# SRA_iris druggability pipeline
# -------------------------------
# Validation case: ERalpha ligand-binding domain (PDB 3ERT), which has the real
# breast-cancer drug tamoxifen (residue code OHT) sitting in its drug pocket.
#
# The plan:
#   1. Download the structure from the PDB.
#   2. Clean it: pull the drug out and save it as an "answer key", and make a
#      bare protein file (no drug, no water) for pocket detection.
#   3. Run fpocket on the bare protein to detect + score pockets.
#   4. (later) Check fpocket's top pocket sits where the drug used to be.
#   5. (later) Visualize in ChimeraX.
#
# Run from anywhere:  bash pipeline.sh
set -euo pipefail

# ---- Config -----------------------------------------------------------------
PDB_ID="3ERT"          # the structure we're using
DRUG_CODE="OHT"        # 3-letter code for tamoxifen inside the file
WORKDIR="$(cd "$(dirname "$0")" && pwd)"   # the folder this script lives in
cd "$WORKDIR"

echo "=============================================================="
echo " SRA_iris pipeline  |  structure: $PDB_ID  |  drug: $DRUG_CODE"
echo " Working folder: $WORKDIR"
echo "=============================================================="

# ---- Step 1: Download -------------------------------------------------------
# Grab the raw structure file from the RCSB PDB. If we already have it, skip.
echo ""
echo "[Step 1] Downloading $PDB_ID.pdb from the PDB ..."
if [[ -f "${PDB_ID}.pdb" ]]; then
    echo "         already have ${PDB_ID}.pdb, skipping download."
else
    curl -sO "https://files.rcsb.org/download/${PDB_ID}.pdb"
    echo "         saved ${PDB_ID}.pdb"
fi

# ---- Step 2: Clean ----------------------------------------------------------
# A PDB file is just text, one atom per line:
#   - lines starting with "ATOM"   = the protein itself
#   - lines starting with "HETATM" = everything else (drug, water, ions)
# The molecule name lives in columns 18-20 (e.g. OHT = drug, HOH = water).
#
# We produce two files:
#   (a) drug_only.pdb     -> just the drug, our "answer key" for later
#   (b) protein_clean.pdb -> only the protein (drug + water removed) for fpocket
echo ""
echo "[Step 2] Cleaning the structure ..."

awk -v drug="$DRUG_CODE" \
    'substr($0,1,6)=="HETATM" && substr($0,18,3)==drug' \
    "${PDB_ID}.pdb" > drug_only.pdb
echo "         wrote drug_only.pdb      (answer key: $(grep -c . drug_only.pdb) drug atoms)"

awk 'substr($0,1,4)=="ATOM"' "${PDB_ID}.pdb" > protein_clean.pdb
echo "         wrote protein_clean.pdb  (bare protein: $(grep -c . protein_clean.pdb) atoms, no drug/water)"

# ---- Step 3: Run fpocket ----------------------------------------------------
# fpocket scans the bare protein, finds every cavity a molecule could bind in,
# and scores each (volume, druggability, polarity, ...).
#
# Command:  fpocket -f protein_clean.pdb
# It creates a folder next to the input: protein_clean_out/  containing
#   - protein_clean_info.txt  -> the numbers for every pocket (the report card)
#   - protein_clean_out.pdb   -> protein + pockets as dummy atoms (for viewing)
#   - pockets/                -> one file per pocket
echo ""
echo "[Step 3] Running fpocket on the bare protein ..."

# Find fpocket whether or not PATH is refreshed in this shell.
if command -v fpocket >/dev/null 2>&1; then
    FPOCKET="fpocket"
elif [[ -x "$HOME/.local/bin/fpocket" ]]; then
    FPOCKET="$HOME/.local/bin/fpocket"
else
    echo "         ERROR: fpocket not found. Install it first." >&2
    exit 1
fi

rm -rf protein_clean_out                        # clear any previous run
"$FPOCKET" -f protein_clean.pdb >/dev/null 2>&1
echo "         fpocket finished. Results in: protein_clean_out/"

NPOCKETS=$(grep -c '^Pocket' protein_clean_out/protein_clean_info.txt)
echo "         pockets found: $NPOCKETS"
echo ""
echo "   --- Top-ranked pocket (Pocket 1) ---"
awk '/^Pocket 1 /{flag=1} flag{print "   "$0} /^Pocket 2 /{exit}' \
    protein_clean_out/protein_clean_info.txt

# ---- Step 4: Does Pocket 1 sit where the drug was? --------------------------
# Both the pocket and the drug are now just clouds of 3D points, so we can
# measure their overlap with actual numbers (no eyeballing needed):
#
#   (a) Center-to-center distance: average position of Pocket 1's marker points
#       vs. average position of the drug atoms. Pockets/drugs are ~10-15 A wide,
#       so if the centers are only a few A apart, they're the same spot.
#   (b) Coverage: of the drug's 29 atoms, how many have a Pocket 1 marker within
#       4 A (i.e. actually sit inside the pocket). High coverage = the drug lives
#       in this pocket.
echo ""
echo "[Step 4] Checking Pocket 1 vs. the drug (answer key) ..."

POCKET1="protein_clean_out/pockets/pocket1_vert.pqr"

awk '
  # Pocket 1 marker points (whitespace fields: x=$6 y=$7 z=$8)
  FILENAME==pf && /^ATOM/ { px[np]=$6; py[np]=$7; pz[np]=$8; np++;
                            sx+=$6; sy+=$7; sz+=$8 }
  # Drug atoms (fixed PDB columns for x/y/z)
  FILENAME==df && /^HETATM/ { x=substr($0,31,8)+0; y=substr($0,39,8)+0; z=substr($0,47,8)+0;
                              dx[nd]=x; dy[nd]=y; dz[nd]=z; nd++;
                              tx+=x; ty+=y; tz+=z }
  END {
    pcx=sx/np; pcy=sy/np; pcz=sz/np           # Pocket 1 center
    dcx=tx/nd; dcy=ty/nd; dcz=tz/nd           # drug center
    dist=sqrt((pcx-dcx)^2+(pcy-dcy)^2+(pcz-dcz)^2)
    cov=0
    for(i=0;i<nd;i++){                        # nearest pocket marker to each drug atom
      best=1e9
      for(j=0;j<np;j++){
        d=sqrt((dx[i]-px[j])^2+(dy[i]-py[j])^2+(dz[i]-pz[j])^2)
        if(d<best) best=d
      }
      if(best<=4.0) cov++
    }
    printf "         Pocket 1 center : (%.1f, %.1f, %.1f)\n", pcx,pcy,pcz
    printf "         Drug center     : (%.1f, %.1f, %.1f)\n", dcx,dcy,dcz
    printf "         Center-to-center distance : %.2f A\n", dist
    printf "         Drug atoms inside Pocket 1 : %d / %d  (%.0f%% within 4 A)\n", cov, nd, 100*cov/nd
    print  ""
    if (dist <= 8.0 && cov >= 0.7*nd)
      print "         VERDICT: Pocket 1 overlaps the drug site -> fpocket found the REAL pocket. PASS."
    else
      print "         VERDICT: Pocket 1 does NOT match the drug site. Investigate."
  }
' pf="$POCKET1" df="drug_only.pdb" "$POCKET1" drug_only.pdb

echo ""
echo "Done through Step 4."

# ---- Step 5: Prepare a Chimera-friendly file for Pocket 1 --------------------
# fpocket wrote Pocket 1's points as a .pqr file with loose (free-format)
# spacing that Chimera can misread. We reformat those points into a clean,
# fixed-column PDB (pocket1.pdb) so Chimera opens them reliably as a blob of
# dummy atoms.
echo ""
echo "[Step 5] Writing pocket1.pdb (Pocket 1's points) for Chimera ..."
awk '/^ATOM/{ n++;
      printf "HETATM%5d  C   STP A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n", \
             n, n, $6, $7, $8 }' \
    "$POCKET1" > pocket1.pdb
echo "         wrote pocket1.pdb ($(grep -c '^HETATM' pocket1.pdb) sphere points)"

echo ""
echo "Done through Step 5. Open view_pocket1.cmd in Chimera to see the result."
