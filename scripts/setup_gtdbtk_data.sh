#!/usr/bin/env bash
set -euo pipefail
# GTDB-Tk DB setup (ACE/UQ primary, ecogenomic.org fallback)
# Modes: latest (full tar) | R226 (split)
# Idempotent: resumes downloads, never deletes or overwrites existing files.
# Now with size verification via Content-Length.

RELEASE="${1:-latest}"   # latest | R226
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEST="${2:-$ROOT/refdata/gtdbtk}"
PKG="gtdbtk_data.tar.gz"

mkdir -p "$DEST"
cd "$DEST"

# ---------------- size helpers ------------------------------------------------
remote_size() {  # remote_size URL -> prints bytes or 0
  local url="$1"
  curl -fsIL -H 'Accept-Encoding: identity' "$url" \
    | awk 'tolower($0) ~ /^content-length:/ {gsub("\r",""); sz=$2} END{print sz+0}'
}

size_close() {  # size_close LOCAL_BYTES REMOTE_BYTES
  local l="$1" r="$2"
  [[ "$r" -le 0 ]] && return 1
  local diff=$(( l>r ? l-r : r-l ))
  # tolerance: max(1%, 1 MiB)
  local tol_pct=$(( r / 100 ))
  local tol_abs=$(( 1024 * 1024 ))
  local tol=$(( tol_pct > tol_abs ? tol_pct : tol_abs ))
  (( diff <= tol ))
}

# ---------------- tolerant download wrappers ---------------------------------
get_file() { # get_file URL [OUT]
  local url="$1" out="${2:-$(basename "$1")}"
  if [[ -f "$out" ]]; then
    local rs=$(remote_size "$url" || echo 0)
    local ls=$(stat -c '%s' "$out" 2>/dev/null || echo 0)
    if size_close "$ls" "$rs"; then
      echo "[gtdbtk] size OK → $out (local $ls / remote $rs)"; return 0
    fi
    echo "[gtdbtk] size mismatch → resume $out (local $ls / remote $rs)"
  fi
  for i in {1..10}; do
    wget -c --retry-connrefused --tries=3 --timeout=60 -O "$out" "$url" && break || true
    sleep $((i*5))
  done
  return 0
}

get_dir() { # recursive, tolerant
  local url="$1"
  wget -c -r -np -nH --cut-dirs=6 -e robots=off -R "index.html*" "$url" || true
}

# ---------------- validators ---------------------------------------------------
is_valid_db() {
  [[ -d auxillary_files ]] || return 1
  compgen -G "release*" >/dev/null || return 1
  [[ -f bac120_taxonomy_r226.tsv ]] || return 1
  [[ -f ar53_taxonomy_r226.tsv    ]] || return 1
  [[ -f bac120_msa_mask_r226.txt  ]] || return 1
  [[ -f ar53_msa_mask_r226.txt    ]] || return 1
  # quick spot checks inside release226
  [[ -f release226/pplacer/gtdb_r226_bac120.refpkg/CONTENTS.json ]] || return 1
  [[ -f release226/markers/pfam/Pfam-A.hmm ]] || return 1
  [[ -f release226/msa/gtdb_r226_bac120.faa ]] || return 1
  [[ -f release226/radii/gtdb_radii.tsv ]] || return 1
  return 0
}

# ---------------- minor helpers ----------------------------------------------
fetch_md5_latest() {
  { curl -fsSL "https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/MD5SUM.txt" \
    || curl -fsSL "https://data.gtdb.ecogenomic.org/releases/latest/MD5SUM.txt"; } \
  | awk '/gtdbtk_package\/full_package\/gtdbtk_data\.tar\.gz/ {print $1; exit}'
}

ensure_aux_only() {
  if compgen -G "release*" >/dev/null && [[ ! -d auxillary_files ]]; then
    echo "[gtdbtk] auxillary_files missing → fetching only aux (resume)"
    get_dir "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/" \
    || get_dir "https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/"
  fi
}

ensure_root_files() {
  [[ -f bac120_taxonomy_r226.tsv ]] || \
    get_file "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_taxonomy_r226.tsv" bac120_taxonomy_r226.tsv \
    || get_file "https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_taxonomy_r226.tsv" bac120_taxonomy_r226.tsv

  [[ -f ar53_taxonomy_r226.tsv ]] || \
    get_file "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_taxonomy_r226.tsv" ar53_taxonomy_r226.tsv \
    || get_file "https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_taxonomy_r226.tsv" ar53_taxonomy_r226.tsv

  [[ -f bac120_msa_mask_r226.txt ]] || \
    get_file "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/bac120_msa_mask_r226.txt" bac120_msa_mask_r226.txt \
    || get_file "https://data.gtdb.ecogenomic.org/releases/release226/226.0/bac120_msa_mask_r226.txt" bac120_msa_mask_r226.txt

  [[ -f ar53_msa_mask_r226.txt ]] || \
    get_file "https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/ar53_msa_mask_r226.txt" ar53_msa_mask_r226.txt \
    || get_file "https://data.gtdb.ecogenomic.org/releases/release226/226.0/ar53_msa_mask_r226.txt" ar53_msa_mask_r226.txt
}

download_full() {
  local primary="https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package"
  local fallback="https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package"
  if [[ -f "$PKG" ]]; then
    # size check vs remote
    local rs=$(remote_size "$primary/$PKG" || remote_size "$fallback/$PKG" || echo 0)
    local ls=$(stat -c '%s' "$PKG" 2>/dev/null || echo 0)
    if size_close "$ls" "$rs"; then
      echo "[gtdbtk] $PKG size OK → skip download."
      return
    fi
    echo "[gtdbtk] $PKG size mismatch → resume download."
  fi
  get_file "$primary/$PKG" "$PKG" || get_file "$fallback/$PKG" "$PKG"
  # optional MD5 note
  local exp act; exp="$(fetch_md5_latest || true)"; act="$(md5sum "$PKG" | awk '{print $1}')"
  [[ -n "${exp:-}" && "$exp" != "$act" ]] && echo "[gtdbtk] MD5 mismatch noted; keeping tar."
}

download_split_r226() {
  local base1="https://data.ace.uq.edu.au/public/gtdb/data/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package"
  local base2="https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/split_package"
  get_dir "$base1/" || get_dir "$base2/"
  if [[ ! -f "$PKG" ]]; then
    shopt -s nullglob
    local parts=(gtdbtk_r226_data.tar.gz.part_*)
    [[ ${#parts[@]} -gt 0 ]] || { echo "[gtdbtk] no split parts found"; return 0; }
    echo "[gtdbtk] assembling parts → $PKG"
    cat gtdbtk_r226_data.tar.gz.part_* > "$PKG" || true
  else
    echo "[gtdbtk] $PKG exists → skip assemble."
  fi
}

# ---------------- main --------------------------------------------------------
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
    if [[ -f "$PKG" ]]; then
      echo "[gtdbtk] extracting (no overwrite)"
      tar -xzf "$PKG" --skip-old-files || echo "[gtdbtk] warn: tar extraction error, will validate instead"
    fi
    ensure_aux_only
    ensure_root_files
    is_valid_db || { echo "[gtdbtk] DB incomplete. Check space or archive."; exit 5; }
  fi
fi

# persist
grep -q 'GTDBTK_DATA_PATH=' "$HOME/.bashrc" 2>/dev/null || echo "export GTDBTK_DATA_PATH=$DEST/release226" >> "$HOME/.bashrc"
export GTDBTK_DATA_PATH="$DEST/release226"
mkdir -p "$ROOT/refdata"
ln -sfn "$DEST/release226" "$ROOT/refdata/gtdbtk"

# smoke test (non-fatal)
if command -v gtdbtk >/dev/null 2>&1; then
  gtdbtk check_install || true
fi

echo "[gtdbtk] ready: $DEST"
