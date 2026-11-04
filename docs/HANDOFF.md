# Handoff: state of branch `claude/adoring-ramanujan-pjh3sn` (PR #5)

Context for the next Claude session working on MeMaMoRe. Written at the end of a
cloud session (4 CPUs, 15 GB RAM, ~30 GB disk) that could only test small parts.

## Already done on this branch (all pushed)

| Commit topic | What |
|---|---|
| `bin/main.sh` | Fixed syntax error (duplicated line in `run_assembly_wrap`) that broke every command; `plasmids_env` path `env/` -> `envs/`; removed duplicate lines |
| Hardcoded paths | `/media/Box/MeMaMoRe` replaced by repo-relative defaults; still overridable via `ROOT`, `REFDATA_BASE`, `GTDBTK_DATA_PATH`, new `VIBRANT_DB` |
| GUI security | Sample names restricted to `[A-Za-z0-9._-]` (leading alnum), module whitelist, cpus/SRR validation, deletes confined to `SRA/`, warning on non-local `--host` |
| Ingest | `--srr ACC` now works (GUI uses it); `env_ingest` installs `sra-tools` + `pigz`; re-running ingest no longer duplicates manifest rows |
| PLASMe DB | Zenodo fallback in `envs/create_env_plasmids.sh` was unreachable under `set -e`; fixed |
| geNomad | Pinned `genomad>=1.12` (1.7 + keras 3.11 crashes); DB <1.9 is re-downloaded; `GENOMAD_SPLITS=N` lowers MMseqs2 RAM (~17.6 GB without); output symlinks + GUI now use geNomad's real file names (`contigs_summary/contigs_plasmid_summary.tsv`) |
| Test data | `scripts/make_mini_testdata.py` |

## Tested and working (cloud session)

Test sample `mini` = reads simulated from phiX174 (5.4 kb), lambda (48.5 kb), plasmid F (99 kb):

```bash
python3 scripts/make_mini_testdata.py            # -> test/mini/ (gitignored)
bin/main.sh env create sra_tools && bin/main.sh env create ingest
bin/main.sh env create qc && bin/main.sh env create assembly
bin/main.sh run ingest --fastq test/mini/manifest.tsv
bin/main.sh run qc --sample mini
bin/main.sh run assembly --sample mini --set global.threads=4
# geNomad step alone (needs genomad_env + refdata/genomad):
GENOMAD_SPLITS=8 micromamba run -n genomad_env bash modules/plasmids/genomad/exec.sh "$PWD" mini 4 0
```

Expected: 3 contigs ~98.0 / 48.5 / 5.4 kb at ~40x; geNomad calls NODE_1 (F) plasmid (~0.99),
NODE_2/NODE_3 (lambda/phiX) virus. Also verified: SRR download (`--srr SRR1770413`), GUI API
shows QC, assembly and geNomad plasmid results.

## Not tested yet (needs a bigger machine)

Binning (MetaBAT2, Binny, COMEBin, MAGScoT, CheckM2, GTDB-Tk), PLASMe (DB 13.6 GB zipped),
MOB-suite, HOTSPOT, ViralVerify, VirSorter2, VIBRANT, CheckV, DeepARG, RGI, AMR++.
Rough needs: GTDB-Tk ~64+ GB RAM and ~140 GB disk; all DBs + envs ~300 GB.

## Known open issues / next steps

1. `bin/main.sh env create all` fails: no umbrella creators `envs/create_env_viruses.sh`
   and `envs/create_env_args.sh` (only per-tool creators exist).
2. Stub/smoke-test mode (placeholder tools writing correctly named dummy outputs) to test
   the whole chain and the GUI in seconds without DBs.
3. GitHub Actions CI: `bash -n` + shellcheck on all `.sh`, `py_compile`/flake8 on `gui/`,
   later the stub run.
4. Missing tools listed in README: VFDB (virulence, e.g. via ABRicate), iPHoP (virus host),
   multi-sample summary.
5. README is outdated (viruses/ARGs are implemented; `env create all` exists).
6. `api/status.jsonl` is tracked although `api/` is in `.gitignore`; pipeline runs modify it,
   so don't commit those changes.

## Notes for running on a local machine

- Needs Linux (Windows: WSL2). micromamba goes to `~/micromamba`
  (`curl -Ls https://micro.mamba.pm/api/micromamba/linux-64/latest | tar -xvj -C ~/micromamba bin/micromamba`).
- Existing setups with an old `genomad_env`: if geNomad fails with
  "Cannot convert ... to a shape", remove the env and rerun `bin/main.sh env create plasmids`.
- Ask before deleting or overwriting anything under `SRA/` or `refdata/` on a machine with real data.
