#!/usr/bin/env bash
set -euo pipefail

module_default_params(){ cat <<'YAML'
modules:
  binning:
    cpus: 8
    min_contig_len: 1500
    metabat2:
      minCv: 0.5
      minContig: 1500
    checkm2:
      completeness_min: 50
      contamination_max: 10
YAML
}

_need(){ command -v "$1" >/dev/null 2>&1 || { echo "Missing tool: $1"; exit 127; }; }
_preflight(){
  _need bash
  _need micromamba
}

run_binning(){
  _preflight

  # ---- inputs ----
  local MANIFEST="$ROOT/SRA/reads/manifest.tsv"
  [[ -s "$MANIFEST" ]] || { echo "[binning] Missing $MANIFEST (run ingest)"; exit 2; }

  : "${GTDBTK_DATA_PATH:=/media/Box/MeMaMoRe/refdata/gtdbtk/release226}"
  [[ -d "$GTDBTK_DATA_PATH" ]] || { echo "[binning] GTDBTK_DATA_PATH invalid: $GTDBTK_DATA_PATH"; exit 2; }

  # threads from merged params helper (prefer python helper; fallback 8)
  CPUS="$(micromamba run -n env_sra_tools python "$ROOT/scripts/merge_params.py" \
           --defaults "$PARAMS_YAML" --get modules.binning.cpus 2>/dev/null || echo 8)"

  # ---- iterate manifest (sample \t fq1 \t fq2?) ----
  tail -n +2 "$MANIFEST" | tr -d '\r' | while IFS=$'\t' read -r SID FQ1 FQ2 || [[ -n "${SID:-}" ]]; do
    [[ -z "${SID:-}" ]] && continue
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue

    # assembly & qc paths
    local ASM_DIR="$ROOT/SRA/assemblies/spades/$SID"
    local QC_DIR="$ROOT/SRA/assemblies/contig_qc/$SID"
    local LOGS="$ROOT/logs"
    mkdir -p "$LOGS"

    # choose contigs (filtered preferred)
    local CONTIGS="$QC_DIR/contigs.filtered.fasta"
    [[ -s "$CONTIGS" ]] || CONTIGS="$ASM_DIR/contigs.len1000.fasta"
    [[ -s "$CONTIGS" ]] || CONTIGS="$ASM_DIR/contigs.fasta"
    [[ -s "$CONTIGS" ]] || { echo "[binning] $SID no contigs found"; continue; }

    # mapping bam (required by several steps)
    local BAM="$ASM_DIR/map/${SID}.sorted.bam"
    [[ -s "$BAM" ]] || { echo "[binning] $SID no BAM (run assembly/contig_qc first)"; continue; }

    echo "[binning] sample=$SID cpus=$CPUS contigs=$(basename "$CONTIGS") bam=$(basename "$BAM")"

    export GTDBTK_DATA_PATH CPUS

    # ---- coverage table ----
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" BAM="$BAM" \
      bash "$ROOT/modules/binning/coverage/exec.sh"

    # ---- MetaBAT2 ----
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" \
      bash "$ROOT/modules/binning/metabat2/exec.sh"

    # ---- Binny (safe cleanup on --force unless resuming) ----
    local BINNY_OUTDIR="$ROOT/SRA/binning/binny/$SID"
    if [[ "${FORCE:-0}" -eq 1 && "${RESUME:-0}" -eq 0 ]]; then
      if [[ -d "$BINNY_OUTDIR" ]]; then
        case "$BINNY_OUTDIR" in
          "$ROOT"/*)
            echo "[binny] --force: removing existing output dir: $BINNY_OUTDIR"
            rm -rf --one-file-system -- "$BINNY_OUTDIR"
            ;;
          *)
            echo "[binny] REFUSE to remove '$BINNY_OUTDIR' (outside ROOT)"; exit 3;;
        esac
      fi
    fi
    mkdir -p "$BINNY_OUTDIR" "$BINNY_OUTDIR/tmp"

    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" BAM="$BAM" FORCE="${FORCE:-0}" \
      bash "$ROOT/modules/binning/binny/exec.sh"

    # ---- COMEBin (best effort) ----
    local GFA="$ASM_DIR/assembly_graph_with_scaffolds.gfa"
    [[ -s "$GFA" ]] || echo "[comebin] $SID missing $GFA (will try without graph)"
    ROOT="$ROOT" SAMPLE="$SID" BAM="$BAM" \
      bash "$ROOT/modules/binning/comebin/exec.sh" || echo "[comebin] $SID skipped or failed"

    # ---- merge candidate bin sets for MAGScoT ----
    local MAGSCOT_DIR="$ROOT/SRA/binning/magscot/$SID"
    mkdir -p "$MAGSCOT_DIR"
    local MAGSCOT_IN="$MAGSCOT_DIR/contigs_to_bin.tsv"
    : > "$MAGSCOT_IN"
    for f in \
      "$ROOT/SRA/binning/metabat2/$SID/contigs_to_bin.with_set.tsv" \
      "$ROOT/SRA/binning/binny/$SID/contigs_to_bin.with_set.tsv" \
      "$ROOT/SRA/binning/comebin/$SID/contigs_to_bin.with_set.tsv"
    do
      [[ -s "$f" ]] && cat "$f" >> "$MAGSCOT_IN"
    done
    [[ -s "$MAGSCOT_IN" ]] || { echo "[magscot] $SID no bin inputs"; continue; }

    # ---- MAGScoT ----
    ROOT="$ROOT" SAMPLE="$SID" MAGS_ROOT="$ROOT/external/MAGScoT" \
      bash "$ROOT/modules/binning/magscot/exec.sh"

    # ---- QC ----
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/qc_checkm2/exec.sh"
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/qc_rrna_trna/exec.sh"

    # ---- Taxonomy ----
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/gtdbtk/exec.sh"
  done
}
