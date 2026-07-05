#!/usr/bin/env bash
#
# FPocket_Script — generalized druggability pipeline
# --------------------------------------------------
# Works on (almost) any protein in the PDB. It:
#   1. Downloads the structure.
#   2. Reports what's inside (chains, possible ligands).
#   3. Cleans it into a bare protein (drops drug, water, ions; and if you pick a
#      chain, drops DNA + extra crystal copies too).
#   4. Runs fpocket to detect + score pockets.
#   5. IF you named a drug, checks whether fpocket's top pocket sits where that
#      drug binds (validation). Otherwise this step is skipped.
#   6. Prepares files + a Chimera script to visualize the result.
#
# USAGE:
#   bash pipeline.sh --pdb ID [--drug CODE] [--chain LETTER]
#
# EXAMPLES:
#   bash pipeline.sh --pdb 3ERT --drug OHT      # ERα: has tamoxifen, single chain (validation runs)
#   bash pipeline.sh --pdb 7VOX --chain A       # FOXA1: no drug, keep chain A (drops DNA; detection only)
#
set -euo pipefail

# ---- Parse command-line options --------------------------------------------
PDB_ID=""; DRUG_CODE=""; CHAIN=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --pdb)   PDB_ID="${2:-}";    shift 2 ;;
    --drug)  DRUG_CODE="${2:-}"; shift 2 ;;
    --chain) CHAIN="${2:-}";     shift 2 ;;
    -h|--help) echo "usage: bash pipeline.sh --pdb ID [--drug CODE] [--chain LETTER]"; exit 0 ;;
    *) echo "unknown option: $1" >&2; exit 1 ;;
  esac
done
if [[ -z "$PDB_ID" ]]; then
  echo "ERROR: --pdb is required.  e.g.  bash pipeline.sh --pdb 3ERT --drug OHT" >&2
  exit 1
fi

# Each protein gets its own folder under runs/ so runs don't overwrite each other.
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="$SCRIPTDIR/runs/$PDB_ID"
mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "=============================================================="
echo " pipeline | PDB: $PDB_ID | drug: ${DRUG_CODE:-<none>} | chain: ${CHAIN:-<all>}"
echo " Output folder: $OUTDIR"
echo "=============================================================="

# ---- Step 1: Download ------------------------------------------------------
echo ""
echo "[Step 1] Downloading $PDB_ID.pdb ..."
if [[ -f "${PDB_ID}.pdb" ]]; then
  echo "         already have ${PDB_ID}.pdb, skipping download."
else
  curl -sf -o "${PDB_ID}.pdb" "https://files.rcsb.org/download/${PDB_ID}.pdb" \
    || { echo "         ERROR: could not download ${PDB_ID}. Check the PDB ID." >&2; exit 1; }
  echo "         saved ${PDB_ID}.pdb"
fi

# ---- Step 2: Peek at what's inside -----------------------------------------
echo ""
echo "[Step 2] What's in this structure:"
echo -n "         protein/nucleic chains present: "
grep '^ATOM' "${PDB_ID}.pdb" | cut -c22 | sort -u | tr '\n' ' '; echo ""
echo "         non-water HETATM molecules (possible drugs / ligands / ions):"
grep '^HETATM' "${PDB_ID}.pdb" | awk '{print substr($0,18,3)}' | sort | uniq -c \
  | awk '$2!="HOH"{printf "           %s  (%s atoms)\n",$2,$1}' | grep . \
  || echo "           (none besides water)"

# ---- Step 3: Clean ---------------------------------------------------------
echo ""
echo "[Step 3] Cleaning ..."

# Protein: keep ATOM records. If a chain was given, keep only that chain
# (this also drops DNA and extra copies, since they live in other chains).
if [[ -n "$CHAIN" ]]; then
  awk -v ch="$CHAIN" 'substr($0,1,4)=="ATOM" && substr($0,22,1)==ch' \
      "${PDB_ID}.pdb" > protein_clean.pdb
  echo "         kept only chain $CHAIN"
else
  awk 'substr($0,1,4)=="ATOM"' "${PDB_ID}.pdb" > protein_clean.pdb
  echo "         kept all protein (ATOM) records"
fi
NATOM=$(grep -c '^ATOM' protein_clean.pdb || true); NATOM=${NATOM:-0}
echo "         wrote protein_clean.pdb ($NATOM atoms; drug/water/ions removed)"

# Safety check: warn if DNA/RNA slipped through.
if grep '^ATOM' protein_clean.pdb \
     | awk '{r=substr($0,18,3); gsub(/ /,"",r); print r}' \
     | grep -qxE 'DA|DT|DG|DC|A|U|G|C'; then
  echo "         ⚠️  WARNING: DNA/RNA residues are still present in the cleaned protein."
  echo "             Re-run with --chain <LETTER> to keep only the protein."
fi

# Drug (answer key): only if the user named a drug code.
HAVE_DRUG=0
if [[ -n "$DRUG_CODE" ]]; then
  if [[ -n "$CHAIN" ]]; then
    awk -v d="$DRUG_CODE" -v ch="$CHAIN" \
        'substr($0,1,6)=="HETATM" && substr($0,18,3)==d && substr($0,22,1)==ch' \
        "${PDB_ID}.pdb" > drug_only.pdb
  else
    awk -v d="$DRUG_CODE" \
        'substr($0,1,6)=="HETATM" && substr($0,18,3)==d' \
        "${PDB_ID}.pdb" > drug_only.pdb
  fi
  NDRUG=$(grep -c '^HETATM' drug_only.pdb || true); NDRUG=${NDRUG:-0}
  if [[ "$NDRUG" -gt 0 ]]; then
    HAVE_DRUG=1
    echo "         wrote drug_only.pdb (answer key: $NDRUG atoms of $DRUG_CODE)"
  else
    echo "         ⚠️  no atoms named '$DRUG_CODE' found — treating as a no-drug run."
    rm -f drug_only.pdb
  fi
else
  echo "         no --drug given → detection-only run (Step 5 will be skipped)"
fi

# ---- Step 4: Run fpocket ---------------------------------------------------
echo ""
echo "[Step 4] Running fpocket ..."
if command -v fpocket >/dev/null 2>&1; then FPOCKET="fpocket"
elif [[ -x "$HOME/.local/bin/fpocket" ]]; then FPOCKET="$HOME/.local/bin/fpocket"
else echo "         ERROR: fpocket not found. Install it first." >&2; exit 1; fi

rm -rf protein_clean_out
"$FPOCKET" -f protein_clean.pdb >/dev/null 2>&1
INFO="protein_clean_out/protein_clean_info.txt"
NPOCKETS=$(grep -c '^Pocket' "$INFO" 2>/dev/null || true); NPOCKETS=${NPOCKETS:-0}
echo "         done — $NPOCKETS pockets found (results in runs/$PDB_ID/protein_clean_out/)"
if [[ "$NPOCKETS" -eq 0 ]]; then
  echo "         no pockets detected; nothing more to do."; exit 0
fi
echo ""
echo "   --- Top-ranked pocket (Pocket 1) ---"
awk '/^Pocket 1 /{f=1} f{print "   "$0} /^Pocket 2 /{exit}' "$INFO"

# ---- Step 5: Validation (only if we have a drug answer key) -----------------
echo ""
POCKET1="protein_clean_out/pockets/pocket1_vert.pqr"
if [[ "$HAVE_DRUG" -eq 1 ]]; then
  echo "[Step 5] Checking Pocket 1 vs. the drug (answer key) ..."
  awk '
    FILENAME==pf && /^ATOM/ { px[np]=$6; py[np]=$7; pz[np]=$8; np++; sx+=$6; sy+=$7; sz+=$8 }
    FILENAME==df && /^HETATM/ { x=substr($0,31,8)+0; y=substr($0,39,8)+0; z=substr($0,47,8)+0;
                                dx[nd]=x; dy[nd]=y; dz[nd]=z; nd++; tx+=x; ty+=y; tz+=z }
    END {
      pcx=sx/np; pcy=sy/np; pcz=sz/np; dcx=tx/nd; dcy=ty/nd; dcz=tz/nd
      dist=sqrt((pcx-dcx)^2+(pcy-dcy)^2+(pcz-dcz)^2)
      cov=0
      for(i=0;i<nd;i++){ best=1e9
        for(j=0;j<np;j++){ e=sqrt((dx[i]-px[j])^2+(dy[i]-py[j])^2+(dz[i]-pz[j])^2); if(e<best)best=e }
        if(best<=4.0)cov++ }
      printf "         Center-to-center distance : %.2f A\n", dist
      printf "         Drug atoms inside Pocket 1 : %d / %d (%.0f%% within 4 A)\n", cov, nd, 100*cov/nd
      if (dist<=8.0 && cov>=0.7*nd) print "         VERDICT: Pocket 1 = the drug site. PASS."
      else print "         VERDICT: Pocket 1 does NOT match the drug site. Investigate."
    }' pf="$POCKET1" df="drug_only.pdb" "$POCKET1" drug_only.pdb
else
  echo "[Step 5] Skipped — no drug provided, so there's no answer key to check against."
  echo "         (Expected for undrugged targets like FOXA1 — you just get the pocket list.)"
fi

# ---- Step 6: Prepare visualization -----------------------------------------
echo ""
echo "[Step 6] Preparing visualization files ..."
awk '/^ATOM/{n++; printf "HETATM%5d  C   STP A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C\n", n,n,$6,$7,$8}' \
    "$POCKET1" > pocket1.pdb
echo "         wrote pocket1.pdb"

# Generate a Chimera script pointing at THIS protein's files (absolute paths).
# Model IDs depend on whether a drug is loaded: #0 protein [, #1 drug], then pocket.
{
  echo "# Auto-generated Chimera script for $PDB_ID"
  echo "open $OUTDIR/protein_clean.pdb"
  [[ "$HAVE_DRUG" -eq 1 ]] && echo "open $OUTDIR/drug_only.pdb"
  echo "open $OUTDIR/pocket1.pdb"
  echo "background solid white"
  echo "~display"
  echo "ribbon #0"
  echo "color gray #0"
  if [[ "$HAVE_DRUG" -eq 1 ]]; then
    echo "display #1"; echo "represent stick #1"; echo "color red #1"     # drug
    echo "display #2"; echo "represent sphere #2"; echo "color orange #2" # pocket
  else
    echo "display #1"; echo "represent sphere #1"; echo "color orange #1" # pocket (no drug)
  fi
  echo "focus"
} > view.cmd
echo "         wrote view.cmd  (open in Chimera:  chimera $OUTDIR/view.cmd)"

TOPDRUG=$(awk '/^Pocket 1 /{f=1} f&&/Druggability/{print $NF; exit}' "$INFO")
echo ""
echo "Done. Summary for $PDB_ID:  pockets=$NPOCKETS  top-druggability=$TOPDRUG"
echo "   files in: runs/$PDB_ID/"
