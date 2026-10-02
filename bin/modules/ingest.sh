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
#          run_ingest --srr SRRxxxxxx      (single accession; uses $SAMPLE)
#   manifest.tsv columns (tab-separated, header required):
#       sample   fq1   fq2
#   fq2 is optional; leave empty for single-end
# -------------------------------------------------------------------
run_ingest(){
  local FASTQ_MANIFEST="" SRR_ACC=""
  local CANONICAL
  CANONICAL="$(yq -r '.modules.ingest.canonical_manifest // "SRA/reads/manifest.tsv"' "$PARAMS_YAML" 2>/dev/null || echo "SRA/reads/manifest.tsv")"

  # Parse flags coming from main (ignore unrelated)
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --fastq) FASTQ_MANIFEST="${2:-}"; shift 2 ;;
      --fastq=*) FASTQ_MANIFEST="${1#--fastq=}"; shift ;;
      --srr) SRR_ACC="${2:-}"; shift 2 ;;
      --srr=*) SRR_ACC="${1#--srr=}"; shift ;;
      --set|--from|--only|--sample) shift 2 || true ;;
      --test|--resume|--force) shift ;;
      *) shift ;;
    esac
  done

  # Single accession (e.g. from the GUI): wrap it into a one-row manifest
  if [[ -n "$SRR_ACC" ]]; then
    [[ -n "${SAMPLE:-}" ]] || _die "--srr requires --sample"
    FASTQ_MANIFEST="$(mktemp)"
    printf "sample\tsrr\n%s\t%s\n" "$SAMPLE" "$SRR_ACC" > "$FASTQ_MANIFEST"
  fi

  [[ -n "$FASTQ_MANIFEST" && -s "$FASTQ_MANIFEST" ]] || _die "Provide --fastq <manifest.tsv> or --srr <accession>"

  # Read and validate header (must have sample; and AT LEAST one of fq1 or srr)
  local hdr
  hdr="$(head -n1 "$FASTQ_MANIFEST" | tr -d '\r')"
  grep -qw "sample" <<<"$hdr" || _die "manifest header missing column: sample"
  if ! (grep -qw "fq1" <<<"$hdr" || grep -qw "srr" <<<"$hdr"); then
    _die "manifest header must contain at least 'fq1' or 'srr'"
  fi

  # Header column indexes
  local col_sample=0 col_fq1=0 col_fq2=0 col_srr=0 col_md5=0
  local i=1
  for col in $hdr; do
    case "$col" in
      sample) col_sample=$i ;;
      fq1)    col_fq1=$i    ;;
      fq2)    col_fq2=$i    ;;
      srr)    col_srr=$i    ;;
      md5)    col_md5=$i    ;;
    esac
    ((i++))
  done

  # Prepare canonical manifest
  local canon_path="$ROOT/$CANONICAL"
  _init_canonical "$canon_path"

  # Iterate rows
  tail -n +2 "$FASTQ_MANIFEST" | tr -d '\r' | \
  awk -F'\t' 'NF{print}' | \
  while IFS=$'\t' read -r -a cols; do
    # bash array is 0-indexed, awk/cut logic was 1-indexed. offset by 1.
    local sample="" fq1="" fq2="" srr="" md5=""
    [[ $col_sample -gt 0 ]] && sample="${cols[$((col_sample-1))]:-}"
    [[ $col_fq1 -gt 0 ]]    && fq1="${cols[$((col_fq1-1))]:-}"
    [[ $col_fq2 -gt 0 ]]    && fq2="${cols[$((col_fq2-1))]:-}"
    [[ $col_srr -gt 0 ]]    && srr="${cols[$((col_srr-1))]:-}"
    [[ $col_md5 -gt 0 ]]    && md5="${cols[$((col_md5-1))]:-}"

    [[ -n "$sample" ]] || continue

    local norm1="" norm2=""

    if [[ -n "$srr" ]]; then
      # --- SRA Download logic ---
      command -v fasterq-dump >/dev/null 2>&1 && command -v pigz >/dev/null 2>&1 \
        || _die "fasterq-dump/pigz missing in env_ingest; recreate it with: bin/main.sh env create ingest"
      echo "[ingest] $sample: downloading SRR $srr ..."
      local outdir="$ROOT/SRA/reads/$sample"
      _mkdirp "$outdir"
      
      # Use fasterq-dump (sra-tools) then compress
      # We assume env_sra_tools is active (which it should be based on main.sh wrapping)
      if [[ ! -f "$outdir/${srr}_1.fastq.gz" && ! -f "$outdir/${sample}_R1.fastq.gz" ]]; then
          fasterq-dump --split-3 --threads "${CPUS:-4}" --outdir "$outdir" "$srr"
          
          # Compress and rename to standardized form
          if [[ -f "$outdir/${srr}_1.fastq" ]]; then
            pigz -p "${CPUS:-4}" "$outdir/${srr}_1.fastq"
            mv "$outdir/${srr}_1.fastq.gz" "$outdir/${sample}_R1.fastq.gz"
            norm1="$outdir/${sample}_R1.fastq.gz"
          fi
          if [[ -f "$outdir/${srr}_2.fastq" ]]; then
            pigz -p "${CPUS:-4}" "$outdir/${srr}_2.fastq"
            mv "$outdir/${srr}_2.fastq.gz" "$outdir/${sample}_R2.fastq.gz"
            norm2="$outdir/${sample}_R2.fastq.gz"
          fi
          # Single end case
          if [[ -f "$outdir/${srr}.fastq" ]]; then
             pigz -p "${CPUS:-4}" "$outdir/${srr}.fastq"
             mv "$outdir/${srr}.fastq.gz" "$outdir/${sample}_R1.fastq.gz"
             norm1="$outdir/${sample}_R1.fastq.gz"
          fi
      else
          echo "[ingest] $sample: SRA outputs already exist, skipping download."
          norm1="$outdir/${sample}_R1.fastq.gz"
          [[ -f "$outdir/${sample}_R2.fastq.gz" ]] && norm2="$outdir/${sample}_R2.fastq.gz"
      fi

    elif [[ -n "$fq1" ]]; then
      # --- Local file logic ---
      [[ -s "$fq1" ]] || _die "$sample: fq1 not found or empty ($fq1)"
      if [[ -n "$fq2" ]]; then
        [[ -s "$fq2" ]] || _die "$sample: fq2 given but file missing/empty ($fq2)"
      fi

      # Optional MD5 Check
      if [[ -n "$md5" ]]; then
        echo "[ingest] $sample: matching MD5 checksum for $fq1 ..."
        # Only check fq1 for now to keep it simple, or split md5 by comma if multiple were supported
        if ! echo "$md5  $fq1" | md5sum -c --status -; then
          _die "$sample: MD5 checksum failed for $fq1"
        fi
      fi

      # Normalize & symlink under SRA/reads/<sample>
      read -r norm1 norm2 <<<"$(_link_normalized "$sample" "$fq1" "$fq2")"
    else
       _die "$sample: Neither fq1 nor srr provided"
    fi

    # Replace any previous row for this sample so re-runs don't duplicate it
    awk -F'\t' -v s="$sample" 'NR==1 || $1!=s' "$canon_path" > "$canon_path.tmp" && mv "$canon_path.tmp" "$canon_path"

    # Append to canonical manifest (fq2/norm2 may be empty)
    if [[ -n "${norm2:-}" ]]; then
      echo -e "${sample}\t${norm1}\t${norm2}" >> "$canon_path"
    else
      echo -e "${sample}\t${norm1}\t"         >> "$canon_path"
    fi

    echo "[ingest] ${sample} -> $(dirname "$norm1")"
  done

  [[ -n "$SRR_ACC" ]] && rm -f "$FASTQ_MANIFEST"
  echo "[ingest] wrote canonical manifest: $canon_path"
}
