#!/usr/bin/env bash
# modules/binning/map_group/exec.sh
# Multi-sample coverage for one assembly: maps the QC'd reads of every member of
# SAMPLE's binning group onto SAMPLE's contigs (bwa-mem2 + samtools). Existing
# BAMs are reused; the sample's own BAM is the contig_qc one. A new group member
# therefore costs only its new mappings.
# Output: SRA/binning/coverage/<s>/group/<member>.bam(.bai), bams.txt (paths,
# one per line, own sample first), members.txt (sorted member list).
# env in: ROOT SAMPLE MEMBERS (space-separated) PARAMS_YAML [CPUS]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"; : "${MEMBERS:?}"
CPUS="${CPUS:-8}"
# shellcheck source=bin/lib/reads.sh
source "$ROOT/bin/lib/reads.sh"
MINLEN="$(cfg global.contig_min_len 1000)"

ASM="$ROOT/SRA/assemblies/spades/$SAMPLE"
REF="$ASM/contigs.len${MINLEN}.fasta"; [[ -s "$REF" ]] || REF="$ASM/contigs.fasta"
OWN="$ASM/map/${SAMPLE}.sorted.bam"
OUT="$ROOT/SRA/binning/coverage/$SAMPLE/group"
LOG="$ROOT/logs/binning_map_group_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"
[[ -s "$REF" && -s "$OWN" ]] || { echo "[map_group] $SAMPLE needs its assembly + contig_qc BAM first"; exit 2; }

if [[ ! -s "${REF}.0123" && ! -s "${REF}.bwt.2bit.64" ]]; then
  micromamba run -n env_mapping_coverm bwa-mem2 index "$REF" >>"$LOG" 2>&1
fi

read -r -a M <<<"$MEMBERS"
: > "$OUT/bams.txt.tmp"
new=0
for R in "$SAMPLE" "${M[@]}"; do
  grep -qxF "$OUT/$R.bam" "$OUT/bams.txt.tmp" && continue      # own sample listed first, once
  B="$OUT/$R.bam"
  if [[ "$R" == "$SAMPLE" ]]; then
    ln -sfn "$OWN" "$B"; ln -sfn "$OWN.bai" "$B.bai"
  else
    rr="$(qc_reads "$R")" || exit 2
    IFS=$'\t' read -r R1 R2 <<<"$rr"
    [[ -s "$R1" ]] || { echo "[map_group] no QC'd reads for group member $R (run qc for it)"; exit 2; }
    if [[ -s "$B" && -s "$B.bai" && "$B" -nt "$REF" && "$B" -nt "$R1" ]]; then
      :   # reuse
    else
      echo "[map_group] $R reads -> $SAMPLE assembly" | tee -a "$LOG"
      if [[ -n "$R2" ]]; then
        micromamba run -n env_mapping_coverm bash -c \
          "bwa-mem2 mem -t $CPUS '$REF' '$R1' '$R2' 2>>'$LOG' | samtools sort -@ $CPUS -o '$B.tmp' -" >>"$LOG" 2>&1
      else
        micromamba run -n env_mapping_coverm bash -c \
          "bwa-mem2 mem -t $CPUS '$REF' '$R1' 2>>'$LOG' | samtools sort -@ $CPUS -o '$B.tmp' -" >>"$LOG" 2>&1
      fi
      mv -f "$B.tmp" "$B"
      micromamba run -n env_mapping_coverm samtools index "$B" >>"$LOG" 2>&1
      new=$((new + 1))
    fi
  fi
  echo "$B" >> "$OUT/bams.txt.tmp"
done
mv -f "$OUT/bams.txt.tmp" "$OUT/bams.txt"
printf '%s\n' "${M[@]}" | sort -u > "$OUT/members.txt"
echo "[map_group] $SAMPLE: $(wc -l < "$OUT/bams.txt") BAMs ($new new, $(( $(wc -l < "$OUT/bams.txt") - 1 - new )) reused)"
