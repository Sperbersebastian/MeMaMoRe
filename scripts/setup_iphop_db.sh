#!/usr/bin/env bash
# iPHoP host database (iPHoP.latest_rw, ~300 GB compressed) into refdata/iphop/.
# The archive is fetched in its ~10 GB chunks; each chunk is md5-checked, streamed
# into tar and deleted, so archive and database never need disk space at the same time.
# Progress/result in refdata/iphop/dl/status (EXTRACTED when done); the viruses module
# skips iPHoP while a download is in progress. A failed run is restarted from scratch.
# Usage: scripts/setup_iphop_db.sh [DB_VERSION] [DEST]
set -uo pipefail
VERSION="${1:-iPHoP.latest_rw}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D="${2:-$ROOT/refdata/iphop}"; DL="$D/dl"; ST="$DL/status"
B="https://portal.nersc.gov/cfs/m342/iphop/db/"
MIN_FREE_GB="${MIN_FREE_GB:-15}"
mkdir -p "$DL" && cd "$DL" || exit 1
curl -sf "${B}${VERSION}.tar.gz.md5" -o db.md5 || { echo "FAILED md5 list ${B}${VERSION}.tar.gz.md5" > "$ST"; exit 1; }
TOTAL=$(awk 'NR==1{print $1}' db.md5)
# line 1: whole archive; further lines (if any): its chunks
if [[ $(awk 'NF>1' db.md5 | wc -l) -gt 1 ]]; then awk 'NR>1 && NF>1{print $1, $2}' db.md5 > parts.md5
else awk 'NR==1{print $1, $2}' db.md5 > parts.md5; fi

# watchdog: stop before the disk is full
( while sleep 60; do
    f=$(df --output=avail -BG "$D" | tail -n1 | tr -dc 0-9)
    if (( f < MIN_FREE_GB )); then echo "FAILED disk full (${f} GB free) $(date)" > "$ST"; kill 0; fi
  done ) &
WD=$!

feed(){
  while read -r md5 name; do
    for a in 1 2 3 4 5 6 7 8 9 10; do
      echo "DOWNLOADING $name (try $a) $(date)" > "$ST"
      curl -sfL -C - --retry 5 -o "$name" "${B}${name}"
      [[ "$(md5sum "$name" 2>/dev/null | cut -d' ' -f1)" == "$md5" ]] && break
      rm -f "$name"; sleep 60
    done
    [[ "$(md5sum "$name" | cut -d' ' -f1)" == "$md5" ]] || { echo "FAILED chunk $name $(date)" > "$ST"; return 1; }
    cat "$name" && rm -f "$name"
  done < parts.md5
}
feed | tee >(md5sum | cut -d' ' -f1 > total.md5) | tar -xz -C "$D"
rc=$?
kill "$WD" 2>/dev/null; sleep 2
if [[ $rc == 0 && "$(cat total.md5)" == "$TOTAL" ]]; then
  echo "EXTRACTED $(date)" > "$ST"; rm -f parts.md5 total.md5
  echo "[iphop-db] done: $(find "$D" -mindepth 1 -maxdepth 1 -type d -name '*_rw')"
else
  grep -q '^FAILED' "$ST" || echo "FAILED rc=$rc $(date)" > "$ST"; exit 1
fi
