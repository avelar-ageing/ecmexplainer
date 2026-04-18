#!/bin/bash
# ==============================================================================
# eleanor_align_script.sh
# RNA-seq alignment pipeline — Salmon quasi-mapping, GRCh38 Ensembl 111
#
# QUICK START
# ───────────
#   1. Edit the USER CONFIGURATION block below (CSV path, output dirs, etc.)
#   2. Prepare your samples CSV (see FORMAT section below)
#   3. Run:  bash eleanor_align_script.sh
#
# The script is safe to restart — any sample whose quant.sf already exists is
# silently skipped.
#
# REQUIREMENTS
# ────────────
#   conda env containing:  sra-tools  salmon  pigz
#   Create it once with:
#     ~/miniconda3/bin/conda create -y -n rnaseq \
#       -c conda-forge -c bioconda --override-channels \
#       sra-tools salmon pigz
#
# CSV FORMAT
# ──────────
# The CSV must have a header row. Required columns:
#
#   label      — unique sample identifier; used as the output directory name
#                (e.g. "GSM1234_treated" or just "sample1")
#   ended      — library layout: "paired" or "single"
#   geo        — study/project ID (used only in progress messages; can be anything)
#
# You also need EITHER:
#   sample     — SRR/ERR/DRR accession; the script downloads the FASTQ from NCBI
#                Multiple accessions for the same label (technical replicates)
#                should appear as separate rows — they are merged into one
#                Salmon run automatically.
#
#   OR:
#   fastq_r1   — absolute path to an already-downloaded R1 (or unpaired) FASTQ
#   fastq_r2   — absolute path to an already-downloaded R2 FASTQ (paired only;
#                leave blank for single-end)
#
# If BOTH `sample` and `fastq_r1` are provided, the local files take priority
# and no download is attempted.
#
# EXAMPLE CSV (SRA download mode)
# ─────────────────────────────────
#   label,ended,geo,sample
#   ctrl_rep1,paired,MyStudy,SRR1234567
#   ctrl_rep2,paired,MyStudy,SRR1234568
#   treated_rep1,paired,MyStudy,SRR1234569
#
# EXAMPLE CSV (local FASTQ mode)
# ─────────────────────────────────
#   label,ended,geo,fastq_r1,fastq_r2
#   ctrl_rep1,paired,MyStudy,/data/ctrl1_R1.fastq.gz,/data/ctrl1_R2.fastq.gz
#   ctrl_rep2,paired,MyStudy,/data/ctrl2_R1.fastq.gz,/data/ctrl2_R2.fastq.gz
#   ctrl_rep3,single,MyStudy,/data/ctrl3.fastq.gz,
#
# EXAMPLE CSV (technical replicates — multiple SRR lanes merged per label)
# ─────────────────────────────────
#   label,ended,geo,sample
#   sampleA,paired,MyStudy,SRR1111111
#   sampleA,paired,MyStudy,SRR1111112
#   sampleA,paired,MyStudy,SRR1111113
#   sampleB,paired,MyStudy,SRR2222221
#
# OUTPUT
# ──────
# For each label, Salmon writes to:
#   <OUTPUT_DIR>/<label>/quant.sf          ← transcript-level TPM + counts
#   <OUTPUT_DIR>/<label>/quant.genes.sf    ← gene-level (if --geneMap used)
#   <OUTPUT_DIR>/<label>/aux_info/         ← mapping stats, etc.
#
# After all samples complete, import with tximport in R:
#   files <- list.files("<OUTPUT_DIR>", "quant.sf", recursive=TRUE, full.names=TRUE)
#   txi   <- tximport(files, type="salmon", tx2gene=tx2gene)
# ==============================================================================

set -uo pipefail

# ==============================================================================
# USER CONFIGURATION — edit these paths before running
# ==============================================================================

# Path to your samples CSV (see format above)
CSV="/path/to/your/samples.csv"

# Where Salmon output directories will be created (one per label)
OUTPUT_DIR="/path/to/your/output/align"

# Temporary directory for downloaded FASTQs (deleted after each sample)
FASTQ_TMP="/path/to/your/fastq_tmp"

# Salmon index directory.
# If it doesn't exist, the script will build it from TRANSCRIPTOME_FA.
# Point this at the same index used for the ECMExplainer training data:
#   /home/ro/Downloads/el_test/hg38_salmon_index
# Or set to a new path to build a fresh index.
SALMON_INDEX="/path/to/hg38_salmon_index"

# Transcriptome FASTA — only needed if the index doesn't exist yet.
# Ensembl 111 GRCh38 cDNA (human):
#   https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# For other organisms/releases change this URL and point SALMON_INDEX to a new path.
TRANSCRIPTOME_FA="/path/to/Homo_sapiens.GRCh38.cdna.all.fa.gz"
TRANSCRIPTOME_URL="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"

# conda environment containing sra-tools, salmon, and pigz
CONDA_ENV="${HOME}/miniconda3/envs/rnaseq"

# CPU threads to use for fasterq-dump, pigz, and salmon
THREADS=6

# ==============================================================================
# END OF USER CONFIGURATION
# ==============================================================================

mkdir -p "$OUTPUT_DIR" "$FASTQ_TMP"

# ── Activate conda env ────────────────────────────────────────────────────────
export PATH="${CONDA_ENV}/bin:${PATH}"

for tool in salmon pigz; do
  if ! command -v "$tool" &>/dev/null; then
    echo "ERROR: $tool not found in ${CONDA_ENV}/bin"
    echo "Create the env with:"
    echo "  ~/miniconda3/bin/conda create -y -n rnaseq -c conda-forge -c bioconda --override-channels sra-tools salmon pigz"
    exit 1
  fi
done

# ── Build Salmon index if needed ──────────────────────────────────────────────
if [ ! -d "$SALMON_INDEX" ]; then
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  Building Salmon index"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  if [ ! -f "$TRANSCRIPTOME_FA" ]; then
    echo "  Downloading transcriptome FASTA..."
    wget -q --show-progress -O "$TRANSCRIPTOME_FA" "$TRANSCRIPTOME_URL"
  fi
  echo "  Indexing (this takes ~10 min)..."
  salmon index \
    --transcripts "$TRANSCRIPTOME_FA" \
    --index       "$SALMON_INDEX" \
    --threads     "$THREADS"
  echo "  ✓ Index built at: $SALMON_INDEX"
else
  echo "✓ Salmon index found: $SALMON_INDEX"
fi

# ── Parse CSV into groups: label → {ended, geo, mode, files[]} ───────────────
# Supports both SRA accession rows and local-FASTQ rows.
# Technical replicates (same label, multiple rows) are merged into one Salmon run.
GROUPS_FILE=$(mktemp)

python3 - "$CSV" "$GROUPS_FILE" <<'PYEOF'
import sys, csv, os
from collections import OrderedDict

csv_path = sys.argv[1]
out_path  = sys.argv[2]

groups = OrderedDict()  # label -> {ended, geo, mode, items:[]}

with open(csv_path, newline='') as fh:
    reader = csv.DictReader(fh)
    cols   = reader.fieldnames or []

    for row in reader:
        label = row.get('label', '').strip()
        if not label:
            continue
        ended = row.get('ended', 'paired').strip()
        geo   = row.get('geo', '').strip()

        # Determine mode: local files take priority over SRA
        fq1 = row.get('fastq_r1', '').strip()
        fq2 = row.get('fastq_r2', '').strip()
        srr = row.get('sample', '').strip()

        if fq1:
            mode = 'local'
            item = f"{fq1}|{fq2}"   # pipe-separated; fq2 may be empty for SE
        elif srr:
            mode = 'sra'
            item = srr
        else:
            print(f"  WARNING: row with label '{label}' has no 'sample' or 'fastq_r1' — skipping.")
            continue

        if label not in groups:
            groups[label] = {'ended': ended, 'geo': geo, 'mode': mode, 'items': []}
        elif groups[label]['mode'] != mode:
            print(f"  WARNING: label '{label}' mixes SRA and local modes — skipping mixed rows.")
            continue

        groups[label]['items'].append(item)

with open(out_path, 'w') as fh:
    for label, info in groups.items():
        items_str = ';;'.join(info['items'])  # ;; separator (paths may contain commas)
        fh.write(f"{label}|{info['ended']}|{info['geo']}|{info['mode']}|{items_str}\n")

print(f"  Parsed {len(groups)} unique labels from {csv_path}")
PYEOF

TOTAL=$(wc -l < "$GROUPS_FILE")
LABEL_NUM=0

while IFS='|' read -r LABEL ENDED GEO MODE ITEMS_STR; do
  LABEL_NUM=$((LABEL_NUM + 1))
  OUT="${OUTPUT_DIR}/${LABEL}"

  # Split items on ;; into an array
  IFS=';;' read -ra ITEMS <<< "$ITEMS_STR"
  N_ITEMS=${#ITEMS[@]}

  echo ""
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
  echo "  [${LABEL_NUM}/${TOTAL}]  ${LABEL}"
  echo "  Study: ${GEO}  |  ${ENDED}-end  |  mode: ${MODE}  |  ${N_ITEMS} file-group(s)"
  echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

  # ── Skip if already done ─────────────────────────────────────────────────
  if [ -f "${OUT}/quant.sf" ]; then
    echo "  ✓ Already aligned — skipping."
    continue
  fi

  # Remove any partial output from a prior failed run
  [ -d "$OUT" ] && rm -rf "$OUT"

  # ── Trap: clean up on error ───────────────────────────────────────────────
  TRAP_ITEMS_STR="$ITEMS_STR"
  trap '
    echo "  ERROR — cleaning up partial output for '"${LABEL}"'";
    rm -rf "'"${OUT}"'";
    IFS=";;" read -ra _TI <<< "'"${TRAP_ITEMS_STR}"'";
    for _item in "${_TI[@]}"; do
      _srr="${_item%%|*}";
      rm -f "'"${FASTQ_TMP}"'/${_srr}"*.fastq "'"${FASTQ_TMP}"'/${_srr}"*.fastq.gz 2>/dev/null || true;
    done
  ' ERR

  # ── Build Salmon ARGS array ───────────────────────────────────────────────
  SALMON_ARGS=()

  if [ "$MODE" = "local" ]; then
    # ── LOCAL FASTQ MODE ─────────────────────────────────────────────────────
    echo "  [1/2] Using local FASTQ files..."
    for ITEM in "${ITEMS[@]}"; do
      FQ1="${ITEM%%|*}"
      FQ2="${ITEM##*|}"

      if [ ! -f "$FQ1" ]; then
        echo "  ✗ File not found: $FQ1 — skipping label."
        trap - ERR
        continue 2
      fi

      if [ "$ENDED" = "paired" ] && [ -n "$FQ2" ] && [ "$FQ1" != "$FQ2" ]; then
        if [ ! -f "$FQ2" ]; then
          echo "  ✗ R2 file not found: $FQ2 — skipping label."
          trap - ERR
          continue 2
        fi
        SALMON_ARGS+=(-1 "$FQ1" -2 "$FQ2")
        echo "    paired: $(basename "$FQ1")  +  $(basename "$FQ2")"
      else
        SALMON_ARGS+=(--unmatedReads "$FQ1")
        echo "    single: $(basename "$FQ1")"
      fi
    done

    echo "  [2/2] Quantifying..."

  else
    # ── SRA DOWNLOAD MODE ────────────────────────────────────────────────────
    echo "  [1/3] Downloading ${N_ITEMS} SRR run(s) from NCBI..."
    for SRR in "${ITEMS[@]}"; do
      echo "    fasterq-dump ${SRR}..."
      rm -f "${FASTQ_TMP}/${SRR}"*.fastq "${FASTQ_TMP}/${SRR}"*.fastq.gz

      SUCCESS=0
      for attempt in 1 2 3 4 5; do
        if fasterq-dump "$SRR" \
            --outdir  "$FASTQ_TMP" \
            --threads "$THREADS" \
            --progress \
            --force \
            --split-3; then
          SUCCESS=1
          break
        else
          echo "    ⚠ Attempt ${attempt}/5 failed for ${SRR} — retrying in 60 s..."
          rm -f "${FASTQ_TMP}/${SRR}"*.fastq "${FASTQ_TMP}/${SRR}"*.fastq.gz
          sleep 60
        fi
      done

      if [ "$SUCCESS" -eq 0 ]; then
        echo "    ✗ fasterq-dump failed 5 times for ${SRR} — skipping label."
        trap - ERR
        continue 2
      fi
    done

    echo "  [2/3] Compressing..."
    for SRR in "${ITEMS[@]}"; do
      pigz -p "$THREADS" "${FASTQ_TMP}/${SRR}"*.fastq
    done

    echo "  [3/3] Quantifying..."

    # Detect layout from first SRR's downloaded files
    FIRST_SRR="${ITEMS[0]}"
    if [ -f "${FASTQ_TMP}/${FIRST_SRR}_2.fastq.gz" ]; then
      echo "    Detected: paired-end"
      for SRR in "${ITEMS[@]}"; do
        SALMON_ARGS+=(-1 "${FASTQ_TMP}/${SRR}_1.fastq.gz")
        SALMON_ARGS+=(-2 "${FASTQ_TMP}/${SRR}_2.fastq.gz")
      done
    elif [ -f "${FASTQ_TMP}/${FIRST_SRR}_1.fastq.gz" ]; then
      echo "    Detected: single-end (_1)"
      for SRR in "${ITEMS[@]}"; do
        SALMON_ARGS+=(--unmatedReads "${FASTQ_TMP}/${SRR}_1.fastq.gz")
      done
    else
      echo "    Detected: single-end (no suffix)"
      for SRR in "${ITEMS[@]}"; do
        SALMON_ARGS+=(--unmatedReads "${FASTQ_TMP}/${SRR}.fastq.gz")
      done
    fi
  fi

  # ── Run Salmon quant (same flags used for the ECMExplainer training data) ──
  mkdir -p "$OUT"

  if [[ " ${SALMON_ARGS[*]} " == *" -1 "* ]] || [[ " ${SALMON_ARGS[*]} " == *" -2 "* ]]; then
    # Paired-end
    salmon quant \
      --index   "$SALMON_INDEX" \
      --libType A \
      "${SALMON_ARGS[@]}" \
      --output  "$OUT" \
      --threads "$THREADS" \
      --seqBias \
      --gcBias  \
      --validateMappings \
      --quiet
  else
    # Single-end (no --gcBias — requires paired reads)
    salmon quant \
      --index        "$SALMON_INDEX" \
      --libType      A \
      "${SALMON_ARGS[@]}" \
      --output       "$OUT" \
      --threads      "$THREADS" \
      --seqBias \
      --validateMappings \
      --quiet
  fi

  # ── Remove downloaded FASTQs immediately to save disk ────────────────────
  if [ "$MODE" = "sra" ]; then
    for SRR in "${ITEMS[@]}"; do
      rm -f "${FASTQ_TMP}/${SRR}"*.fastq.gz "${FASTQ_TMP}/${SRR}.fastq.gz"
    done
  fi

  echo "  ✓ Done → ${OUT}/quant.sf"
  trap - ERR

done < "$GROUPS_FILE"

rm -f "$GROUPS_FILE"

echo ""
echo "✅ All ${TOTAL} labels processed."
echo "   Output: ${OUTPUT_DIR}/"
echo ""
echo "Next step — import in R with tximport:"
echo "  files <- list.files('${OUTPUT_DIR}', 'quant.sf', recursive=TRUE, full.names=TRUE)"
echo "  names(files) <- basename(dirname(files))"
echo "  txi <- tximport::tximport(files, type='salmon', tx2gene=tx2gene)"
