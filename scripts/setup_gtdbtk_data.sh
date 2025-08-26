#!/usr/bin/env bash
set -euo pipefail
# GTDB-Tk DB setup (ACE/UQ primary, ecogenomic.org fallback)
# Modes: latest (full tar) | R226 (split)
# Idempotent: resumes downloads, never deletes or overwrites existing files.

RELEASE="${1:-latest}"   # latest | R226
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${2:-$ROOT/refdata/gtdbtk}"
PKG="gtdbtk_data.tar.gz"

mkdir -p "$DEST"
cd "$DEST"

# ---- validators --------------------------------------------------------------
is_valid_db() {
  [[ -d auxillary_files ]] || return 1
  compgen -G "release*" >/dev/null || return 1
  # R226: required files at DB root
  [[ -f bac120_taxonomy_r226.tsv ]] || return 1
  [[ -f ar53_taxonomy_r226.tsv    ]] || return 1
  [[ -f bac120_msa_mask_r226.txt  ]] || return 1
  [[ -f ar53_msa_mask_r226.txt    ]] || return 1
  return 0
}

# ---- helpers ----------------------------------------------------------------
fetch_md5_latest() {
  { curl -fsSL "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM.txt" \
    || curl -fsSL "https://data.gtdb.ecogenomic.org/releases/latest/MD5SUM.txt"; } \
  | awk '/gtdbtk_package\/full_package\/gtdbtk_data\.tar\.gz/ {print $1; exit}'
}

ensure_aux_only() {
  if compgen -G "release*" >/dev/null && [[ ! -d auxillary_files ]]; then
    echo "[gtdbtk] auxillary_files missing → fetching only aux (resume)"
    wget -c -r -np -nH --cut-dirs=6 -e robots=off -R "index.html*" \
      "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/" \
      || wget -c -r -np -nH --cut-dirs=6 -e robots=off -R "index.html*" \
         "https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/"
  fi
}

ensure_root_files() {
  # Always ACE/UQ first, ecogenomic.org fallback
  [[ -f bac120_taxonomy_r226.tsv ]] || wget -c -O bac120_taxonomy_r226.tsv \
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_taxonomy_r226.tsv" \
    || wget -c -O bac120_taxonomy_r226.tsv \
       "https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv"

  [[ -f ar53_taxonomy_r226.tsv ]] || wget -c -O ar53_taxonomy_r226.tsv \
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_taxonomy_r226.tsv" \
    || wget -c -O ar53_taxonomy_r226.tsv \
       "https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_taxonomy_r226.tsv"

  [[ -f bac120_msa_mask_r226.txt ]] || wget -c -O bac120_msa_mask_r226.txt \
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_msa_mask_r226.txt" \
    || wget -c -O bac120_msa_mask_r226.txt \
       "https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_msa_mask_r226.txt"

  [[ -f ar53_msa_mask_r226.txt ]] || wget -c -O ar53_msa_mask_r226.txt \
    "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_msa_mask_r226.txt" \
    || wget -c -O ar53_msa_mask_r226.txt \
       "https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_msa_mask_r226.txt"
}

download_full() {
  local primary="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package"
  local fallback="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package"
  if [[ -f "$PKG" ]]; then
    echo "[gtdbtk] $PKG exists → skip download."
    return
  fi
  wget -c -t 3 --timeout=30 "$primary/$PKG" || wget -c "$fallback/$PKG"
  # Soft MD5 note (no deletion)
  local exp act; exp="$(fetch_md5_latest || true)"; act="$(md5sum "$PKG" | awk '{print $1}')"
  [[ -n "${exp:-}" && "$exp" != "$act" ]] && echo "[gtdbtk] MD5 mismatch noted; keeping tar."
}

download_split_r226() {
  local base1="https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package"
  local base2="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package"
  wget -c -r -nd -np -e robots=off -A 'gtdbtk_r226_data.tar.gz.part_*' "$base1" \
    || wget -c -r -nd -np -e robots=off -A 'gtdbtk_r226_data.tar.gz.part_*' "$base2"
  if [[ ! -f "$PKG" ]]; then
    shopt -s nullglob
    local parts=(gtdbtk_r226_data.tar.gz.part_*)
    [[ ${#parts[@]} -gt 0 ]] || { echo "[gtdbtk] no split parts found"; exit 3; }
    echo "[gtdbtk] assembling parts → $PKG"
    cat gtdbtk_r226_data.tar.gz.part_* > "$PKG"
  else
    echo "[gtdbtk] $PKG exists → skip assemble."
  fi
}

# ---- main -------------------------------------------------------------------
if is_valid_db; then
  echo "[gtdbtk] valid DB present → skip download/extract."
else
  ensure_aux_only
  ensure_root_files
  if ! is_valid_db; then
    case "$RELEASE" in
      latest|LATEST) download_full ;;
      R226|226)      download_split_r226 ;;
      *) echo "Unsupported release: $RELEASE"; exit 2 ;;
    esac
    # extract missing files only
    if [[ -f "$PKG" ]]; then
      echo "[gtdbtk] extracting (no overwrite)"
      tar -xzf "$PKG" --skip-old-files
    fi
    ensure_aux_only
    ensure_root_files
    is_valid_db || { echo "[gtdbtk] DB incomplete. Check space or archive."; exit 5; }
  fi
fi

# Persist env + stable repo link
grep -q 'GTDBTK_DATA_PATH=' "$HOME/.bashrc" 2>/dev/null || echo "export GTDBTK_DATA_PATH=$DEST" >> "$HOME/.bashrc"
export GTDBTK_DATA_PATH="$DEST"
mkdir -p "$ROOT/refdata"
ln -sfn "$DEST" "$ROOT/refdata/gtdbtk"

# Smoke test (non-fatal)
if command -v gtdbtk >/dev/null 2>&1; then
  gtdbtk check_install || true
fi

echo "[gtdbtk] ready: $DEST"
