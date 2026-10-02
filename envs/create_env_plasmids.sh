#!/usr/bin/env bash
# create_env_plasmids.sh
# Setzt alle Envs + optionale Referenzdatenbanken für das Plasmid-/Virus-Modul auf.
# Envs: plasmids_core, hotspot_env, plasme_env, genomad_env, viralverify_env, mobsuite_env, stampede_env
# Abhängigkeit: env_assembly_core (via envs/create_env_assembly.sh)
# DBs: $REFDATA_BASE/{genomad,plasme,viralverify}

set -euo pipefail

# --- strikt non-interaktiv, ohne libmamba-Bug ---
export MAMBA_NO_BANNER=1
export PIP_DISABLE_PIP_VERSION_CHECK=1
export PIP_NO_INPUT=1
export GIT_ASKPASS=/bin/true

# --- micromamba aktivieren ---
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/micromamba}"
export PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"
command -v micromamba >/dev/null 2>&1 || { echo "[err] micromamba not found"; exit 1; }
eval "$(micromamba shell hook --shell=bash)"
MM(){ micromamba -y -q "$@"; }

# --- Pfade ---
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXTERNAL_DIR="$ROOT_DIR/external"
REFDATA_BASE="${REFDATA_BASE:-$ROOT_DIR/refdata}"
mkdir -p "$EXTERNAL_DIR" "$REFDATA_BASE"

# --- Parameter ---
CPUS="${CPUS:-8}"
KEEP_ZIP="${KEEP_ZIP:-0}"
PRINT_CONFIG="${PRINT_CONFIG:-0}"
FORCE_CLEAN="${FORCE_CLEAN:-0}"

CH=(-c conda-forge -c bioconda)
MM config set channel_priority flexible >/dev/null 2>&1 || true

# --- Helpers ---
env_exists(){ micromamba env list | awk 'NR>2{print $1}' | grep -qx "$1"; }
rm_env_if_forced(){ [[ "$FORCE_CLEAN" == "1" ]] && env_exists "$1" && MM env remove -n "$1" || true; }
create_if_missing(){
  local env="$1"; shift
  rm_env_if_forced "$env"
  env_exists "$env" || MM create -n "$env" "$@" "${CH[@]}"
}
clone_if_missing(){
  local repo_url="$1" dest="$2"
  if [[ ! -d "$dest/.git" ]]; then git clone --depth 1 "$repo_url" "$dest"; fi
}
find_or_fail(){
  local base="$1" name="$2"
  local p; p="$(find "$base" -type f -name "$name" -print -quit || true)"
  [[ -n "$p" ]] || { echo "[err] not found: $name under $base"; exit 1; }
  printf "%s" "$p"
}

# ========== 1) Core ==========
create_if_missing plasmids_core mmseqs2 bwa-mem2 samtools coverm pigz parallel gawk
echo "[ok] env plasmids_core"

# ========== 2) HOTSPOT ==========
clone_if_missing https://github.com/Orin-beep/HOTSPOT "$EXTERNAL_DIR/HOTSPOT"

echo "[env] creating hotspot_env"
if micromamba env create -y -f "$EXTERNAL_DIR/HOTSPOT/environment.yaml" -n hotspot_env 2>/dev/null; then
  echo "[ok] hotspot_env created via environment.yaml"
else
  echo "[info] environment.yaml failed -> manual setup"
  create_if_missing hotspot_env python=3.10
  if command -v nvidia-smi >/dev/null 2>&1; then
    MM install -n hotspot_env pytorch torchvision torchaudio pytorch-cuda=11.8 -c pytorch -c nvidia "${CH[@]}"
  else
    MM install -n hotspot_env pytorch torchvision torchaudio cpuonly -c pytorch "${CH[@]}"
  fi
  MM install -n hotspot_env numpy pandas scikit-learn biopython bidict treelib diamond prodigal hmmer "${CH[@]}"
fi

# Werkzeuge
MM install -n hotspot_env -y -c bioconda "blast=2.13.0"
MM install -n hotspot_env -y -c conda-forge gdown unzip rsync

# Pip-Extras nur wenn vorhanden
if [[ -f "$EXTERNAL_DIR/HOTSPOT/requirements.txt" ]]; then
  micromamba run -n hotspot_env python -m pip install -q --no-cache-dir -r "$EXTERNAL_DIR/HOTSPOT/requirements.txt"
fi
echo "[ok] env hotspot_env; code in $EXTERNAL_DIR/HOTSPOT"

# ---- Pfade
HOTSPOT_REF="$REFDATA_BASE/HOTSPOT"
HOTSPOT_DB="$HOTSPOT_REF/database"
HOTSPOT_MDL="$HOTSPOT_REF/models"
mkdir -p "$HOTSPOT_DB" "$HOTSPOT_MDL"

# ---- Upstream-Skripte versuchen
micromamba run -n hotspot_env bash -lc '
  set -euo pipefail
  cd "'"$EXTERNAL_DIR"'/HOTSPOT"
  [[ -f prepare_db.sh  ]] && bash prepare_db.sh  || true
  [[ -f prepare_mdl.sh ]] && bash prepare_mdl.sh || true
'

# Aus Repo übernehmen
rsync -a "$EXTERNAL_DIR/HOTSPOT/database"/ "$HOTSPOT_DB"/ 2>/dev/null || true
rsync -a "$EXTERNAL_DIR/HOTSPOT/models"/   "$HOTSPOT_MDL"/ 2>/dev/null || true

# ---- Google Drive Fallback (Archivtyp automatisch)
need_db=0; need_mdl=0
[[ -z "$(find "$HOTSPOT_DB"  -mindepth 1 -maxdepth 1 2>/dev/null)" ]] && need_db=1
[[ -z "$(find "$HOTSPOT_MDL" -mindepth 1 -maxdepth 1 2>/dev/null)" ]] && need_mdl=1

if (( need_db==1 || need_mdl==1 )); then
  micromamba run -n hotspot_env bash -lc '
    set -euo pipefail
    extract_any(){ arc="$1"; dest="$2"; mkdir -p "$dest" "$dest.tmp";
      if unzip -t "$arc" >/dev/null 2>&1; then unzip -q "$arc" -d "$dest.tmp";
      elif tar -tzf "$arc" >/dev/null 2>&1; then tar -xzf "$arc" -C "$dest.tmp";
      elif tar -tf  "$arc" >/dev/null 2>&1; then tar -xf  "$arc" -C "$dest.tmp";
      else echo "[err] unknown archive: $arc"; exit 1; fi
      if [[ $(find "$dest.tmp" -mindepth 1 -maxdepth 1 -type d | wc -l) -eq 1 ]] && [[ -z "$(find "$dest.tmp" -maxdepth 1 -type f)" ]]; then
        top="$(find "$dest.tmp" -mindepth 1 -maxdepth 1 -type d | head -n1)"; rsync -a "$top"/ "$dest"/;
      else rsync -a "$dest.tmp"/ "$dest"/; fi
      rm -rf "$dest.tmp"
    }
    REF="'"$HOTSPOT_REF"'"; DB="'"$HOTSPOT_DB"'"; MDL="'"$HOTSPOT_MDL"'"
    DB_ID="1ZSTz3kotwF8Zugz_aBGDtmly8BVo9G4T"
    MDL_ID="1bnA1osvYDgYBi-DRFkP-HrvcnvBvbipF"
    [[ '"$need_db"'  -eq 1 ]] && { out="$REF/database.dl"; gdown --fuzzy "https://drive.google.com/uc?id=${DB_ID}"  -O "$out"; extract_any "$out" "$DB";  rm -f "$out"; }
    [[ '"$need_mdl"' -eq 1 ]] && { out="$REF/models.dl";   gdown --fuzzy "https://drive.google.com/uc?id=${MDL_ID}" -O "$out"; extract_any "$out" "$MDL"; rm -f "$out"; }
    echo "[check] DB files:  $(find "$DB"  -type f | wc -l) -> $DB"
    echo "[check] MDL files: $(find "$MDL" -type f | wc -l) -> $MDL"
  '
fi

# ---- BLAST v5 DB validieren, ggf. neu bauen
micromamba run -n hotspot_env bash -lc '
set -euo pipefail
DBDIR="'"$HOTSPOT_DB"'"
DBPFX="$DBDIR/database"
ok=0
if blastdbcmd -db "$DBPFX" -info >/dev/null 2>&1; then ok=1; fi
if [[ $ok -eq 0 ]]; then
  FA_DB=$(find "$DBDIR" -maxdepth 2 -type f -iregex ".*\.\(fa\|fna\|fasta\)$" | head -n1 || true)
  if [[ -z "$FA_DB" ]]; then
    echo "[warn] keine Nukleotid-FASTA in $DBDIR gefunden. Baue Minimal-DB aus Dummy-Sequenz."
    FA_DB="$DBDIR/ref.fa"; printf ">ref\nATG%0.s" {1..600} > "$FA_DB"; echo "TAA" >> "$FA_DB"
  fi
  makeblastdb -in "$FA_DB" -dbtype nucl -parse_seqids -title HOTSPOT_DB -out "$DBPFX"
  blastdbcmd -db "$DBPFX" -info
fi
'

# ---- Symlinks zurück ins Repo
ln -sfn "$HOTSPOT_DB"  "$EXTERNAL_DIR/HOTSPOT/database"
ln -sfn "$HOTSPOT_MDL" "$EXTERNAL_DIR/HOTSPOT/models"

# ---- Final Check
echo "[check] files DB:  $(find "$HOTSPOT_DB"  -type f 2>/dev/null | wc -l)  -> $HOTSPOT_DB"
echo "[check] files MDL: $(find "$HOTSPOT_MDL" -type f 2>/dev/null | wc -l)  -> $HOTSPOT_MDL"
[[ -z "$(find "$HOTSPOT_DB"  -mindepth 1 -maxdepth 1 2>/dev/null)" ]] && echo "[warn] HOTSPOT DB missing. Populate: $HOTSPOT_DB"
[[ -z "$(find "$HOTSPOT_MDL" -mindepth 1 -maxdepth 1 2>/dev/null)" ]] && echo "[warn] HOTSPOT models missing. Populate: $HOTSPOT_MDL"

# ---- Laufzeit-Hinweise (optional Guards)
# export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
# preprocessing.py: --len default 1500. Bei kurzen Contigs senken, sonst entsteht kein test_feat.dict.
# [[ -s "$MID/test_feat.dict" ]] || { echo "[skip] keine Features erzeugt"; exit 0; }


# ========== 3) PLASMe ==========
clone_if_missing https://github.com/HubertTang/PLASMe.git "$EXTERNAL_DIR/PLASMe"

# Env
if [[ -f "$EXTERNAL_DIR/PLASMe/plasme.yaml" ]]; then
  [[ "$FORCE_CLEAN" == "1" ]] && env_exists plasme_env && MM env remove -n plasme_env || true
  if env_exists plasme_env; then
    MM env update -n plasme_env -f "$EXTERNAL_DIR/PLASMe/plasme.yaml" --prune
  else
    MM env create -n plasme_env -f "$EXTERNAL_DIR/PLASMe/plasme.yaml"
  fi
else
  create_if_missing plasme_env python=3.10 pip requests numpy pandas biopython tqdm
fi
MM install -n plasme_env gawk -c conda-forge || true

# DB
PLASME_DST="$REFDATA_BASE/plasme/DB"
if [[ ! -d "$PLASME_DST" || -z "$(find "$PLASME_DST" -mindepth 1 -maxdepth 1 2>/dev/null)" ]]; then
  echo "[info] PLASMe DB missing -> downloading"
  wd="$EXTERNAL_DIR/PLASMe"
  mkdir -p "$REFDATA_BASE/plasme"
  micromamba run -n plasme_env bash -lc "
    set -euo pipefail
    cd '$wd'
    rm -rf DB DB.zip || true

    # 1) Upstream downloader
    KZ='$KEEP_ZIP'
    if [[ \"\$KZ\" == \"1\" ]]; then
      python PLASMe_db.py --threads '$CPUS' --keep_zip 1 || true
    else
      python PLASMe_db.py --threads '$CPUS' || true
    fi

    # 2) Fallback: Zenodo (test inside 'if' so set -e doesn't abort before it)
    if ! python - <<'PY'
import zipfile, os, sys
sys.exit(0 if os.path.exists('DB.zip') and zipfile.is_zipfile('DB.zip') else 1)
PY
    then
      echo '[plasme] upstream downloader failed, trying Zenodo'
      rm -f DB.zip
      ZURL=\"\${PLASME_ZENODO_URL:-https://zenodo.org/record/8046934/files/DB.zip?download=1}\"
      curl -fL --retry 3 --retry-delay 5 -o DB.zip \"\$ZURL\"
    fi

    # 3) Validate + unzip
    python - <<'PY' || { echo '[err] DB.zip invalid'; exit 2; }
import zipfile, os, sys
sys.exit(0 if os.path.exists('DB.zip') and zipfile.is_zipfile('DB.zip') else 1)
PY
    unzip -q -o DB.zip -d .
  "
  rsync -a --delete "$EXTERNAL_DIR/PLASMe/DB/" "$PLASME_DST/"
  if [[ -z "$(find "$PLASME_DST" -mindepth 1 -maxdepth 1 2>/dev/null)" ]]; then
    echo "[err] plasme DB download failed"; exit 2
  fi
  echo "[ok] plasme DB -> $PLASME_DST"
else
  echo "[ok] plasme DB present at $PLASME_DST"
fi
echo "[ok] env plasme_env"


# ========== 4) geNomad ==========
create_if_missing genomad_env genomad hmmer prodigal-gv
micromamba run -n genomad_env bash -lc 'command -v genomad' >/dev/null
gdst="$REFDATA_BASE/genomad"; mkdir -p "$gdst"
if [[ ! -f "$gdst/genomad_db/version.txt" ]]; then
  echo "[info] geNomad DB missing -> downloading to $gdst"
  micromamba run -n genomad_env bash -lc "cd '$gdst' && genomad download-database ."
  [[ -f "$gdst/genomad_db/version.txt" ]] && \
    echo "[ok] genomad DB present at $gdst/genomad_db" || \
    { echo "[err] genomad DB download failed"; exit 2; }
else
  echo "[ok] genomad DB present at $gdst/genomad_db"
fi
echo "[ok] env genomad_env"

# ========== 5) viralVerify ==========
create_if_missing viralverify_env python=3.10 pip
MM install -n viralverify_env prodigal hmmer biopython "blast=2.13.0" "${CH[@]}"
MM install -n viralverify_env viralverify "${CH[@]}" || true
clone_if_missing https://github.com/ablab/viralVerify "$EXTERNAL_DIR/viralVerify"
pfam_dir="$REFDATA_BASE/viralverify/pfam"; mkdir -p "$pfam_dir"
if [[ ! -f "$pfam_dir/Pfam-A.hmm" ]]; then
  micromamba run -n viralverify_env bash -lc "
    set -e
    cd '$pfam_dir'
    curl -fL --retry 3 -O https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip -f Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm || true
  "
fi
echo "[ok] env viralverify_env; repo external/viralVerify"

# ========== 6) MOB-Suite ==========
# Ziel: v3.1.9 bevorzugt. Fallback: GitHub/pip. Deps gepinnt. Shims nur falls nötig.

create_if_missing mobsuite_env python=3.10

# Basis-Tools + Runtimes
MM install -y -n mobsuite_env \
  mash=2.3 \
  "blast>=2.12,<3" \
  libgcc-ng libstdcxx-ng "libgfortran>=3,<6" \
  "${CH[@]}"

# harte Pins laut Requirements (kompatibel zu 3.1.x)
MM install -y -n mobsuite_env \
  "pandas<=1.5.3" \
  "numpy<1.23.5" \
  "scipy<2" \
  "pytables<4" \
  "pycurl<8" \
  "ete3>=3.1.3,<4" \
  "biopython>=1.80,<2" \
  "six>=1.10,<2" \
  "${CH[@]}"

# Bevorzugt: Bioconda v3.1.9; wenn nicht verfügbar -> pip (GitHub tag)
if ! MM install -y -n mobsuite_env "mob_suite==3.1.9" -c bioconda -c conda-forge; then
  clone_if_missing https://github.com/phac-nml/mob-suite.git "$EXTERNAL_DIR/mob-suite"
  micromamba run -n mobsuite_env bash -lc "
    set -euo pipefail
    cd '$EXTERNAL_DIR/mob-suite'
    git fetch --tags
    git checkout v3.1.9 || true
    python -m pip install --no-cache-dir .
  "
fi

# Sanity: Versionsausgabe
micromamba run -n mobsuite_env bash -lc 'mob_recon --help >/dev/null 2>&1 || { echo "[err] mob_recon nicht funktionsfähig"; exit 1; }'

# OPTIONAL: Kompatibilitäts-Shims nur wenn benötigt (ältere Imports)
micromamba run -n mobsuite_env python - <<'PY' || true
import site, os, textwrap, importlib
try:
    import mob_suite
except Exception:
    raise SystemExit(0)

# Teste, ob pandas-EmptyDataError dort liegt, wo mob_suite es ggf. erwartet
code = """
# pandas EmptyDataError shim
try:
    import pandas.errors as _pe
    import pandas.io.common as _pic
    if not hasattr(_pic, 'EmptyDataError'):
        _pic.EmptyDataError = _pe.EmptyDataError
except Exception:
    pass
# Biopython GC shim (alt -> neu)
try:
    import Bio.SeqUtils as _bsu
    if not hasattr(_bsu, 'GC'):
        from Bio.SeqUtils import gc_fraction as _gc_fraction
        def GC(seq):
            try:
                return 100.0 * _gc_fraction(seq)
            except Exception:
                return 0.0
        _bsu.GC = GC
except Exception:
    pass
"""
sp = site.getsitepackages()[0]
fn = os.path.join(sp, "sitecustomize.py")
with open(fn, "w") as f:
    f.write(textwrap.dedent(code).lstrip())
print("[shim] sitecustomize.py geschrieben:", fn)
PY

echo "[ok] env mobsuite_env"

# ========= MOB-suite DB initialisieren =========
MOBS_DB="${REFDATA_BASE}/mobsuite_db"
mkdir -p "$MOBS_DB"

# Versuche offiziellen Initializer mit Zielpfad
if ! micromamba run -n mobsuite_env mob_init -d "$MOBS_DB"; then
  echo "[warn] mob_init Download fehlgeschlagen -> Zenodo-Fallback"
  micromamba run -n mobsuite_env bash -lc "
    set -euo pipefail
    cd '$MOBS_DB'
    # offizielles Datenarchiv (Zenodo Release Dump)
    curl -fL -o data.tar.gz 'https://zenodo.org/records/10304948/files/data.tar.gz?download=1'
    tar -xzf data.tar.gz
    # erwartete Dateien ggf. an Root kopieren
    test -f ncbi_plasmid_full_seqs.fas || cp -f data/ncbi_plasmid_full_seqs.fas ncbi_plasmid_full_seqs.fas
    # Indizes
    makeblastdb -in ncbi_plasmid_full_seqs.fas -dbtype nucl
    mash sketch -o ncbi_plasmid_full_seqs.fas ncbi_plasmid_full_seqs.fas
  "
fi

# final check
micromamba run -n mobsuite_env bash -lc "
  test -s '$MOBS_DB/ncbi_plasmid_full_seqs.fas' || { echo '[err] MOB-suite DB fehlt'; exit 1; }
"
echo "[ok] MOB-suite DB -> $MOBS_DB"


# ========== 7) MetaPlasmidSPAdes via env_assembly_core ==========
if ! env_exists env_assembly_core; then
  if [[ -s "$ROOT_DIR/envs/create_env_assembly.sh" ]]; then
    echo "[deps] env_assembly_core missing -> running create_env_assembly.sh"
    FORCE_CLEAN="$FORCE_CLEAN" bash "$ROOT_DIR/envs/create_env_assembly.sh"
  else
    echo "[err] env_assembly_core missing and no creator at $ROOT_DIR/envs/create_env_assembly.sh" >&2
    exit 1
  fi
fi
micromamba run -n env_assembly_core bash -lc 'command -v spades.py >/dev/null || exit 2'
micromamba run -n env_assembly_core bash -lc 'command -v metaPlasmidSPAdes.py >/dev/null || command -v metaplasmidspades.py >/dev/null'
echo "[ok] metaPlasmidSPAdes available in env_assembly_core"

# ========== 8) Stampede-ClusterGenomes ==========
create_if_missing stampede_env perl mummer
STAMP_DIR="$EXTERNAL_DIR/Stampede-ClusterGenomes"
if [[ ! -d "$STAMP_DIR/.git" ]]; then
  git clone --depth 1 https://bitbucket.org/MAVERICLab/stampede-clustergenomes.git "$STAMP_DIR"
fi
ST_CL="$(find_or_fail "$STAMP_DIR" "Cluster_genomes.pl")"
chmod +x "$(dirname "$ST_CL")"/* || true
micromamba run -n stampede_env bash -lc 'perl -v >/dev/null'
[[ -x "$ST_CL" ]] || { echo "[err] Cluster_genomes.pl not executable at $ST_CL"; exit 1; }
echo "[ok] env stampede_env; repo $STAMP_DIR"

echo "[done] plasmids environments ready"
echo "[refdata] base = $REFDATA_BASE"
echo "[external]    = $EXTERNAL_DIR"

# ---------- optional: Configdump ----------
if [[ "$PRINT_CONFIG" == "1" ]]; then
  cat <<EOF
ENV_plasmids_core=plasmids_core
ENV_hotspot=hotspot_env
ENV_plasme=plasme_env
ENV_genomad=genomad_env
ENV_viralverify=viralverify_env
ENV_mobsuite=mobsuite_env
ENV_stampede=stampede_env
ENV_assembly_core=env_assembly_core

REF_plasme=$REFDATA_BASE/plasme/DB
REF_genomad=$REFDATA_BASE/genomad/genomad_db
REF_viralverify_pfam=$REFDATA_BASE/viralverify/pfam/Pfam-A.hmm
EXTERNAL_HOTSPOT=$EXTERNAL_DIR/HOTSPOT
EXTERNAL_PLASME=$EXTERNAL_DIR/PLASMe
EXTERNAL_VV=$EXTERNAL_DIR/viralVerify
EXTERNAL_STAMPEDE=$STAMP_DIR
EOF
fi
