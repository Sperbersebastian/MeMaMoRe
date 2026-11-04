#!/usr/bin/env bash
# modules/binning/semibin2/exec.sh
# SemiBin2 (self-supervised) as 4th binner; its bins go into MAGScoT. Best
# effort like COMEBin: bin/modules/binning.sh continues when it fails.
# Uses the central gene calls (genes.faa) instead of its own ORF finding.
# BAMs: BAMS (space-separated; several = multi-sample coverage, B16) or the
# sample's own contig_qc BAM. Output: SRA/binning/semibin2/<s>/output_bins/*.fa
# and contigs_to_bin.with_set.tsv (contig, bin, semibin2).
# env in: ROOT SAMPLE [CPUS] [FORCE] [PARAMS_YAML] [BAMS]
set -euo pipefail

: "${ROOT:?}"; : "${SAMPLE:?}"
CPUS="${CPUS:-8}"; FORCE="${FORCE:-0}"
OUT="$ROOT/SRA/binning/semibin2/$SAMPLE"
LOG="$ROOT/logs/binning_semibin2_${SAMPLE}.log"
mkdir -p "$OUT" "$ROOT/logs"
# shellcheck source=bin/lib/cfg.sh
source "$ROOT/bin/lib/cfg.sh"
case "$(cfg modules.binning.semibin2 true)" in
  false|False|0|no) echo "[semibin2] disabled (modules.binning.semibin2: false)"; exit 0 ;;
esac
ENVIRONMENT="$(cfg modules.binning.semibin2_environment "")"   # empty = self-supervised

CONTIGS="$ROOT/SRA/assemblies/spades/$SAMPLE/contigs.fasta"
GENES="$ROOT/SRA/assemblies/genes/$SAMPLE/genes.faa"
read -r -a BAM_LIST <<<"${BAMS:-$ROOT/SRA/assemblies/spades/$SAMPLE/map/${SAMPLE}.sorted.bam}"
[[ -s "$CONTIGS" ]] || { echo "[semibin2] no contigs for $SAMPLE"; exit 2; }
for b in "${BAM_LIST[@]}"; do [[ -s "$b" ]] || { echo "[semibin2] missing BAM $b"; exit 2; }; done

MAP="$OUT/contigs_to_bin.with_set.tsv"
if [[ -s "$MAP" && "$FORCE" != "1" && "$MAP" -nt "${BAM_LIST[0]}" ]]; then
  echo "[semibin2] cached -> $MAP"; exit 0
fi

micromamba env list --json | grep -q '/semibin_env"' || bash "$ROOT/envs/create_env_semibin.sh"

rm -rf "$OUT/output_bins"
cmd=(SemiBin2 single_easy_bin -i "$CONTIGS" -b "${BAM_LIST[@]}" -o "$OUT" -t "$CPUS"
     --compression none --engine cpu)
if [[ -n "$ENVIRONMENT" ]]; then cmd+=(--environment "$ENVIRONMENT"); else cmd+=(--self-supervised); fi
[[ -s "$GENES" ]] && cmd+=(--prodigal-output-faa "$GENES")
echo "[semibin2] ${#BAM_LIST[@]} BAM(s): ${cmd[*]}" >>"$LOG"
micromamba run -n semibin_env "${cmd[@]}" >>"$LOG" 2>&1

: > "$MAP.tmp"
shopt -s nullglob
for fa in "$OUT"/output_bins/*.fa; do
  b="$(basename "${fa%.fa}")"
  awk -v b="$b" '/^>/{id=substr($1,2); print id "\t" b "\tsemibin2"}' "$fa" >> "$MAP.tmp"
done
shopt -u nullglob
mv -f "$MAP.tmp" "$MAP"
echo "[semibin2] $(cut -f2 "$MAP" | sort -u | wc -l) bins -> $OUT/output_bins"
