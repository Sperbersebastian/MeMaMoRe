#!/usr/bin/env bash
# Shared lookup of the QC'd reads every downstream step should use.
#   qc_reads <sample>   -> prints "R1<TAB>R2" (R2 empty for single-end)
# Default: fastp output (SRA/qc/fastp/<s>/<s>_trimmed_{1,2}.fastq.gz).
# With modules.qc.host_removal = human: Hostile output (SRA/qc/hostile/<s>/
# <s>_clean_{1,2}.fastq.gz); missing host-cleaned reads are an error, so
# unfiltered reads are never used by accident.
# Needs ROOT; reads PARAMS_YAML through bin/lib/cfg.sh.

# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"

qc_reads(){
  local s="$1" mode r1 r2
  mode="$(cfg modules.qc.host_removal none)"
  case "$mode" in
    none)
      r1="$ROOT/SRA/qc/fastp/$s/${s}_trimmed_1.fastq.gz"
      r2="$ROOT/SRA/qc/fastp/$s/${s}_trimmed_2.fastq.gz" ;;
    human)
      r1="$ROOT/SRA/qc/hostile/$s/${s}_clean_1.fastq.gz"
      r2="$ROOT/SRA/qc/hostile/$s/${s}_clean_2.fastq.gz"
      [[ -s "$r1" ]] || { echo "[reads] host_removal=human but no host-cleaned reads for $s ($r1); run: bin/main.sh run qc --sample $s" >&2; return 2; } ;;
    *) echo "[reads] invalid modules.qc.host_removal '$mode' (none|human)" >&2; return 2 ;;
  esac
  [[ -s "$r2" ]] || r2=""
  printf '%s\t%s\n' "$r1" "$r2"
}
