#!/usr/bin/env python3
"""Turn a simulation config (config/sim.yaml format) into plain tables for modules/sim.

    sim_plan.py <config.yaml> <outdir>

Writes to <outdir>:
  settings.tsv    key/value: set_name, n_reads, read_len, seed, model, cpus
  references.tsv  accession, type (host|plasmid|phage), host (accession or empty)
  weights.tsv     sample, accession, weight (fraction of the sample's read pairs)
  provenance.tsv  type, accession, set, host  (the truth table used by tests/e2e/evaluate.py)

Two config formats are accepted.

Legacy (config/sim.yaml, config/testsets/*.yaml):
  hosts: {hq: ACC, lq: ACC}   plasmids: [ACC, ...]   phages: [ACC, ...]
  abundances: {ACC: weight, ...}
  -> one sample named after set_name; no host links.

Extended:
  references:
    - {acc: NC_000913.3, type: host, abundance: 0.30}
    - {acc: NC_013353.1, type: plasmid, host: NC_000913.3, abundance: 0.007}
  samples:                      # optional, explicit per-sample weights
    - {name: a, abundances: {NC_000913.3: 0.2, ...}}
  profiles: {n: 6, sigma: 1.0, copy_sigma: 0.2}   # optional, generated samples

With `profiles`, sample i gets every independent reference's base abundance multiplied by
exp(N(0, sigma)). A plasmid or phage with a host follows its host's factor times
exp(N(0, copy_sigma)), so element and host abundances co-vary across samples as they do in
real communities. Each sample is rescaled to the total weight of the base abundances.
Sample names are <set_name>_s1 ... <set_name>_sN. The random generator is seeded with `seed`.
"""
import math
import random
import sys
from pathlib import Path

import yaml

TYPES = ("host", "plasmid", "phage")


def _refs_legacy(cfg):
    hosts = cfg.get("hosts") or {}
    if isinstance(hosts, dict):
        hosts = [v for v in hosts.values() if v]
    refs = [{"acc": a, "type": "host"} for a in hosts]
    refs += [{"acc": a, "type": "plasmid"} for a in (cfg.get("plasmids") or [])]
    refs += [{"acc": a, "type": "phage"} for a in (cfg.get("phages") or [])]
    abund = cfg.get("abundances") or {}
    for r in refs:
        r["abundance"] = float(abund.get(r["acc"], 0.0))
    return refs


def load_references(cfg):
    if "references" in cfg:
        refs = []
        for r in cfg["references"]:
            refs.append({"acc": str(r["acc"]), "type": r["type"], "host": r.get("host") or "",
                         "abundance": float(r.get("abundance", 0.0))})
    else:
        refs = _refs_legacy(cfg)
        for r in refs:
            r["host"] = ""
    seen = set()
    for r in refs:
        if r["type"] not in TYPES:
            raise SystemExit(f"[sim_plan] {r['acc']}: type must be one of {TYPES}, got {r['type']!r}")
        if r["acc"] in seen:
            raise SystemExit(f"[sim_plan] duplicate accession {r['acc']}")
        seen.add(r["acc"])
    hosts = {r["acc"] for r in refs if r["type"] == "host"}
    for r in refs:
        if r["host"] and r["host"] not in hosts:
            raise SystemExit(f"[sim_plan] {r['acc']}: host {r['host']} is not listed as type host")
    return refs


def make_samples(cfg, refs, set_name, seed):
    base = {r["acc"]: r["abundance"] for r in refs}
    if cfg.get("samples"):
        out = []
        for s in cfg["samples"]:
            w = {acc: float(s.get("abundances", {}).get(acc, 0.0)) for acc in base}
            out.append((str(s["name"]), w))
        return out
    prof = cfg.get("profiles")
    if not prof:
        return [(set_name, base)]
    n = int(prof.get("n", 1))
    sigma = float(prof.get("sigma", 1.0))
    copy_sigma = float(prof.get("copy_sigma", 0.2))
    rng = random.Random(seed)
    total = sum(base.values())
    out = []
    for i in range(1, n + 1):
        factor = {r["acc"]: math.exp(rng.gauss(0.0, sigma)) for r in refs if not r["host"]}
        for r in refs:
            if r["host"]:
                factor[r["acc"]] = factor[r["host"]] * math.exp(rng.gauss(0.0, copy_sigma))
        w = {acc: base[acc] * factor[acc] for acc in base}
        s = sum(w.values())
        if s > 0:
            w = {acc: v * total / s for acc, v in w.items()}
        out.append((f"{set_name}_s{i}", w))
    return out


def main(argv):
    if len(argv) != 3:
        raise SystemExit(__doc__)
    cfg = yaml.safe_load(Path(argv[1]).read_text()) or {}
    out = Path(argv[2])
    out.mkdir(parents=True, exist_ok=True)
    set_name = str(cfg.get("set_name") or "").strip()
    if not set_name:
        raise SystemExit("[sim_plan] set_name missing")
    seed = int(cfg.get("seed", 42))
    refs = load_references(cfg)
    samples = make_samples(cfg, refs, set_name, seed)

    for k in ("n_reads", "read_len"):
        if k not in cfg:
            raise SystemExit(f"[sim_plan] {k} missing")
    with open(out / "settings.tsv", "w") as fh:
        fh.write(f"set_name\t{set_name}\nn_reads\t{int(cfg['n_reads'])}\n"
                 f"read_len\t{int(cfg['read_len'])}\nseed\t{seed}\n"
                 f"model\t{cfg.get('model', 'HiSeq')}\ncpus\t{int(cfg.get('cpus', 8))}\n")
    with open(out / "references.tsv", "w") as fh:
        fh.write("accession\ttype\thost\n")
        for r in refs:
            fh.write(f"{r['acc']}\t{r['type']}\t{r['host']}\n")
    with open(out / "weights.tsv", "w") as fh:
        fh.write("sample\taccession\tweight\n")
        for name, w in samples:
            for acc in (r["acc"] for r in refs):
                fh.write(f"{name}\t{acc}\t{w[acc]:.6f}\n")
    with open(out / "provenance.tsv", "w") as fh:
        fh.write("type\taccession\tset\thost\n")
        for r in refs:
            fh.write(f"{r['type']}\t{r['acc']}\t{set_name}\t{r['host']}\n")
    print(f"[sim_plan] {set_name}: {len(refs)} references, {len(samples)} sample(s) -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
