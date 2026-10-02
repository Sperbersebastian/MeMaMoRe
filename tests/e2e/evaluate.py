#!/usr/bin/env python3
"""Score MeMaMoRe plasmid predictions against the ground truth of a simulated test set.

Inputs come from tests/e2e/run_e2e.sh:
  --refs        combined reference FASTA used for read simulation (refdata/testsets/<SET>/refs/combined.fna)
  --provenance  truth table written by modules/sim (type, accession, set)
  --contigs-paf minimap2 alignment of the assembly contigs against --refs
  --pred-paf    minimap2 alignment of the plasmid union (headers "<tool>|<name>") against --refs

Every contig / predicted sequence is labelled with the reference it aligns to best
(most aligned bases); sequences whose best reference covers less than --min-qcov of
their length are "unassigned" and count as false positives.

Metrics per tool (plus "union" over all tools):
  precision_bp       plasmid bases / all predicted bases
  precision_n        plasmid sequences / all predicted sequences
  recall_bp_assembly reference-plasmid bases covered by predictions / covered by assembly contigs
  recall_bp_ref      reference-plasmid bases covered by predictions / total plasmid length
  f1_bp              harmonic mean of precision_bp and recall_bp_assembly
"""
import argparse
import csv
import json
import sys
from collections import defaultdict


def read_fasta_lengths(path):
    lengths, name, n = {}, None, 0
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line)
    if name is not None:
        lengths[name] = n
    return lengths


def read_provenance(path):
    types = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            types[row["accession"]] = row["type"]
    return types


def read_paf(path):
    """Yield (qname, qlen, tname, tstart, tend, aligned_bases) per alignment line."""
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                continue
            qname, qlen = f[0], int(f[1])
            qstart, qend = int(f[2]), int(f[3])
            tname, tstart, tend = f[5], int(f[7]), int(f[8])
            yield qname, qlen, tname, tstart, tend, qend - qstart


def assign(paf_path, min_len, min_qcov):
    """Return {query: (qlen, best_target_or_None, [(target, start, end), ...])}."""
    per_q = {}
    for qname, qlen, tname, tstart, tend, alen in read_paf(paf_path):
        if qlen < min_len:
            continue
        q = per_q.setdefault(qname, {"qlen": qlen, "bases": defaultdict(int), "iv": []})
        q["bases"][tname] += alen
        q["iv"].append((tname, tstart, tend))
    out = {}
    for qname, q in per_q.items():
        best, bases = max(q["bases"].items(), key=lambda kv: kv[1])
        target = best if bases >= min_qcov * q["qlen"] else None
        out[qname] = (q["qlen"], target, [iv for iv in q["iv"] if iv[0] == target])
    return out


def merged_length(intervals):
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e:
            if cur_e is not None:
                total += cur_e - cur_s
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s
    return total


def covered_plasmid_bases(assigned, types):
    by_ref = defaultdict(list)
    for _qlen, target, ivs in assigned.values():
        if target is not None and types.get(target) == "plasmid":
            for t, s, e in ivs:
                by_ref[t].append((s, e))
    return sum(merged_length(v) for v in by_ref.values())


def safe_div(a, b):
    return round(a / b, 4) if b else None


def score(tool, assigned, types, asm_cov, ref_total):
    n = len(assigned)
    bp = sum(qlen for qlen, _t, _iv in assigned.values())
    tp = [(qlen, t) for qlen, t, _iv in assigned.values() if t is not None and types.get(t) == "plasmid"]
    covered = covered_plasmid_bases(assigned, types)
    p_bp = safe_div(sum(q for q, _ in tp), bp)
    r_asm = safe_div(covered, asm_cov)
    f1 = round(2 * p_bp * r_asm / (p_bp + r_asm), 4) if p_bp and r_asm else None
    return {
        "tool": tool,
        "n_predicted": n,
        "bp_predicted": bp,
        "n_plasmid": len(tp),
        "precision_n": safe_div(len(tp), n),
        "precision_bp": p_bp,
        "recall_bp_assembly": r_asm,
        "recall_bp_ref": safe_div(covered, ref_total),
        "f1_bp": f1,
    }


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--refs", required=True)
    ap.add_argument("--provenance", required=True)
    ap.add_argument("--contigs-paf", required=True)
    ap.add_argument("--pred-paf", required=True)
    ap.add_argument("--min-len", type=int, default=1000, help="ignore contigs shorter than this (bp)")
    ap.add_argument("--min-qcov", type=float, default=0.5, help="aligned fraction needed to assign a label")
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-json")
    args = ap.parse_args(argv)

    types = read_provenance(args.provenance)
    ref_len = read_fasta_lengths(args.refs)
    ref_total = sum(l for acc, l in ref_len.items() if types.get(acc) == "plasmid")

    contigs = assign(args.contigs_paf, args.min_len, args.min_qcov)
    asm_cov = covered_plasmid_bases(contigs, types)

    preds = assign(args.pred_paf, args.min_len, args.min_qcov)
    by_tool = defaultdict(dict)
    for qname, v in preds.items():
        tool = qname.split("|", 1)[0] if "|" in qname else "unknown"
        by_tool[tool][qname] = v

    rows = [score(t, by_tool[t], types, asm_cov, ref_total) for t in sorted(by_tool)]
    rows.append(score("union", preds, types, asm_cov, ref_total))

    cols = list(rows[0].keys()) if rows else []
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if v is None else v) for k, v in r.items()})

    summary = {
        "plasmid_ref_bp": ref_total,
        "plasmid_bp_in_assembly": asm_cov,
        "assembly_recall_bp_ref": safe_div(asm_cov, ref_total),
        "tools": rows,
    }
    if args.out_json:
        with open(args.out_json, "w") as fh:
            json.dump(summary, fh, indent=2)
    json.dump(summary, sys.stdout, indent=2)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())
