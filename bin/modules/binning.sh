#!/usr/bin/env bash
set -euo pipefail

module_default_params(){ cat <<'YAML'
modules:
  binning:
    min_contig_len: 1500
    metabat2:
      minCv: 0.5
      minContig: 1500
    checkm2:
      completeness_min: 50
      contamination_max: 10
YAML
}

run_binning(){
  local SAMPLES_TSV="$ROOT/config/samples.tsv"
  [[ -s "$SAMPLES_TSV" ]] || { echo "Missing $SAMPLES_TSV"; exit 2; }

  tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r SID COUNTRY MGMT R1 R2 FASTA; do
    [[ -n "${SAMPLE:-}" && "$SID" != "$SAMPLE" ]] && continue

    ASM_DIR="$ROOT/SRA/assemblies/spades/$SID"
    QC_DIR="$ROOT/SRA/assemblies/contig_qc/$SID"
    BINDIR="$ROOT/SRA/binning"; LOGS="$ROOT/logs"
    mkdir -p "$BINDIR"/{coverage,metabat2,binny,comebin,magscot,checkm2,barrnap,trnascan,gtdbtk}/"$SID" "$LOGS"

    CONTIGS="$QC_DIR/contigs.filtered.fasta"
    [[ -s "$CONTIGS" ]] || CONTIGS="$ASM_DIR/contigs.len1000.fasta"
    [[ -s "$CONTIGS" ]] || { echo "[binning] $SID no contigs"; continue; }

    BAM="$ASM_DIR/map/${SID}.sorted.bam"
    [[ -s "$BAM" ]] || { echo "[binning] $SID no BAM (run assembly/contig_qc first)"; continue; }

    # coverage
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" BAM="$BAM" \
      bash "$ROOT/modules/binning/coverage/exec.sh"

    # metabat2
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" \
      bash "$ROOT/modules/binning/metabat2/exec.sh"

    # binny
    ROOT="$ROOT" PARAMS_YAML="$PARAMS_YAML" SAMPLE="$SID" CONTIGS="$CONTIGS" BAM="$BAM" \
      bash "$ROOT/modules/binning/binny/exec.sh"

    # COMEBin
    # uses SPAdes graph + contigs.fasta for IDs, BAM for depth (computes depth if missing)
    GFA="$ASM_DIR/assembly_graph_with_scaffolds.gfa"
    [[ -s "$GFA" ]] || { echo "[comebin] $SID missing $GFA"; }
    ROOT="$ROOT" SAMPLE="$SID" BAM="$BAM" \
      bash "$ROOT/modules/binning/comebin/exec.sh"

    # combine for MAGScoT
    MAGSCOT_DIR="$ROOT/SRA/binning/magscot/$SID"
    mkdir -p "$MAGSCOT_DIR"
    MAGSCOT_IN="$MAGSCOT_DIR/contigs_to_bin.tsv"
    : > "$MAGSCOT_IN"
    for f in \
      "$ROOT/SRA/binning/metabat2/$SID/contigs_to_bin.with_set.tsv" \
      "$ROOT/SRA/binning/binny/$SID/contigs_to_bin.with_set.tsv" \
      "$ROOT/SRA/binning/comebin/$SID/contigs_to_bin.with_set.tsv"
    do
      [[ -s "$f" ]] && cat "$f" >> "$MAGSCOT_IN"
    done
    [[ -s "$MAGSCOT_IN" ]] || { echo "[magscot] $SID no bin inputs"; continue; }

    # magscot
    ROOT="$ROOT" SAMPLE="$SID" MAGS_ROOT="$ROOT/tools/MAGScoT" \
      bash "$ROOT/modules/binning/magscot/exec.sh"

    # qc
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/qc_checkm2/exec.sh"
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/qc_rrna_trna/exec.sh"

    # taxonomy
    ROOT="$ROOT" SAMPLE="$SID" bash "$ROOT/modules/binning/gtdbtk/exec.sh"
  done
}
