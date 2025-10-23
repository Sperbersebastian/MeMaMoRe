#!/usr/bin/env bash
set -euo pipefail
[[ "${DEBUG:-0}" == "1" ]] && set -x

# -------------------------------------------------------------------
# Default params (kept tiny on purpose; no metadata)
# -------------------------------------------------------------------
module_default_params(){ cat <<'YAML'
modules:
  ingest:
    # write a canonical manifest with normalized paths
    canonical_manifest: "SRA/reads/manifest.tsv"
YAML
}

# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------
_die(){ echo "[ingest] $*" >&2; exit 2; }

# Ensure directory exists
_mkdirp(){ mkdir -p "$1"; }

# Write header if canonical manifest missing
_init_canonical(){
  local canon="$1"
  local canon_dir
  canon_dir="$(dirname "$canon")"
  _mkdirp "$canon_dir"
  [[ -s "$canon" ]] || echo -e "sample\tfq1\tfq2" > "$canon"
}

# Normalize into SRA/reads/<sample>/sample_R{1,2}.fastq[.gz] (symlinks)
_link_normalized(){
  local sample="$1" fq1="$2" fq2="${3:-}"
  local outdir="$ROOT/SRA/reads/$sample"
  _mkdirp "$outdir"

  # always create R1
  ln -sfn "$fq1" "$outdir/${sample}_R1.fastq.gz" 2>/dev/null || ln -sfn "$fq1" "$outdir/${sample}_R1.fastq"
  # optional R2
  if [[ -n "${fq2:-}" ]]; then
    ln -sfn "$fq2" "$outdir/${sample}_R2.fastq.gz" 2>/dev/null || ln -sfn "$fq2" "$outdir/${sample}_R2.fastq"
    echo "$outdir/${sample}_R1.fastq.gz" "$outdir/${sample}_R2.fastq.gz" >/dev/null 2>&1 || true
  fi

  # return paths we actually wrote (prefer .gz names)
  local p1="$outdir/${sample}_R1.fastq.gz"; [[ -e "$p1" ]] || p1="$outdir/${sample}_R1.fastq"
  local p2=""
  if [[ -n "${fq2:-}" ]]; then
    p2="$outdir/${sample}_R2.fastq.gz"; [[ -e "$p2" ]] || p2="$outdir/${sample}_R2.fastq"
  fi
  printf "%s\t%s\n" "$p1" "$p2"
}

# -------------------------------------------------------------------
# Entry point
#   Usage: run_ingest --fastq manifest.tsv
#   manifest.tsv columns (tab-separated, header required):
#       sample   fq1   fq2
#   fq2 is optional; leave empty for single-end
# -------------------------------------------------------------------
run_ingest(){
  local FASTQ_MANIFEST=""
  local CANONICAL
  CANONICAL="$(yq -r '.modules.ingest.canonical_manifest // "SRA/reads/manifest.tsv"' "$PARAMS_YAML" 2>/dev/null || echo "SRA/reads/manifest.tsv")"

  # Parse flags coming from main (ignore unrelated)
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fastq) FASTQ_MANIFEST="${2:-}"; shift 2 ;;
      --fastq=*) FASTQ_MANIFEST="${1#--fastq=}"; shift ;;
      --set|--from|--only|--sample) shift 2 || true ;;
      --test|--resume|--force) shift ;;
      *) shift ;;
    esac
  done

  [[ -n "$FASTQ_MANIFEST" && -s "$FASTQ_MANIFEST" ]] || _die "Provide --fastq <manifest.tsv>"

  # Read and validate header
  local hdr
  hdr="$(head -n1 "$FASTQ_MANIFEST" | tr -d '\r')"
  grep -qw "sample" <<<"$hdr" || _die "manifest header missing column: sample"
  grep -qw "fq1"    <<<"$hdr" || _die "manifest header missing column: fq1"
  # fq2 optional

  # Prepare canonical manifest
  local canon_path="$ROOT/$CANONICAL"
  _init_canonical "$canon_path"

  # Iterate rows
  tail -n +2 "$FASTQ_MANIFEST" | tr -d '\r' | \
  awk -F'\t' 'NF{print}' | \
  while IFS=$'\t' read -r sample fq1 fq2 rest; do
    [[ -n "${sample:-}" ]] || continue

    [[ -s "${fq1:-}" ]] || _die "$sample: fq1 not found or empty ($fq1)"
    if [[ -n "${fq2:-}" ]]; then
      [[ -s "$fq2" ]] || _die "$sample: fq2 given but file missing/empty ($fq2)"
    fi

    # Normalize & symlink under SRA/reads/<sample>
    read -r norm1 norm2 <<<"$(_link_normalized "$sample" "$fq1" "${fq2:-}")"

    # Append to canonical manifest (fq2 may be empty)
    if [[ -n "${norm2:-}" ]]; then
      echo -e "${sample}\t${norm1}\t${norm2}" >> "$canon_path"
    else
      echo -e "${sample}\t${norm1}\t"         >> "$canon_path"
    fi

    echo "[ingest] ${sample} -> $(dirname "$norm1")"
  done

  echo "[ingest] wrote canonical manifest: $canon_path"
}
