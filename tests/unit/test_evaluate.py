"""Unit tests for tests/e2e/evaluate.py on a tiny hand-built truth set."""
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("evaluate", HERE.parent / "e2e" / "evaluate.py")
evaluate = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evaluate)


def paf(*rows):
    # qname qlen qstart qend strand tname tlen tstart tend matches alen mapq
    return "".join(
        f"{q}\t{ql}\t{qs}\t{qe}\t+\t{t}\t{tl}\t{ts}\t{te}\t{qe - qs}\t{qe - qs}\t60\n"
        for q, ql, qs, qe, t, tl, ts, te in rows
    )


def build(tmp_path):
    (tmp_path / "refs.fna").write_text(">CHR\n" + "A" * 10000 + "\n>PL1\n" + "C" * 4000 + "\n>PL2\n" + "G" * 2000 + "\n")
    (tmp_path / "prov.tsv").write_text("type\taccession\tset\nhost\tCHR\tt\nplasmid\tPL1\tt\nplasmid\tPL2\tt\n")
    # assembly recovers PL1 fully (two contigs), PL2 not at all, plus one chromosome contig
    (tmp_path / "contigs.paf").write_text(paf(
        ("c1", 2500, 0, 2500, "PL1", 4000, 0, 2500),
        ("c2", 1500, 0, 1500, "PL1", 4000, 2500, 4000),
        ("c3", 5000, 0, 5000, "CHR", 10000, 0, 5000),
    ))
    # genomad: c1 correct + c3 false positive; plasme: c2 correct; one short seq ignored
    (tmp_path / "pred.paf").write_text(paf(
        ("genomad|c1", 2500, 0, 2500, "PL1", 4000, 0, 2500),
        ("genomad|c3", 5000, 0, 5000, "CHR", 10000, 0, 5000),
        ("plasme|c2", 1500, 0, 1500, "PL1", 4000, 2500, 4000),
        ("plasme|tiny", 500, 0, 500, "PL2", 2000, 0, 500),
    ))


def run(tmp_path):
    build(tmp_path)
    out = tmp_path / "m.json"
    evaluate.main([
        "--refs", str(tmp_path / "refs.fna"), "--provenance", str(tmp_path / "prov.tsv"),
        "--contigs-paf", str(tmp_path / "contigs.paf"), "--pred-paf", str(tmp_path / "pred.paf"),
        "--out-tsv", str(tmp_path / "m.tsv"), "--out-json", str(out),
    ])
    data = json.loads(out.read_text())
    return data, {r["tool"]: r for r in data["tools"]}


def test_assembly_level_truth(tmp_path):
    data, _ = run(tmp_path)
    assert data["plasmid_ref_bp"] == 6000
    assert data["plasmid_bp_in_assembly"] == 4000
    assert data["assembly_recall_bp_ref"] == round(4000 / 6000, 4)


def test_per_tool_scores(tmp_path):
    _, t = run(tmp_path)
    assert t["genomad"]["precision_bp"] == round(2500 / 7500, 4)
    assert t["genomad"]["recall_bp_assembly"] == round(2500 / 4000, 4)
    assert t["plasme"]["n_predicted"] == 1  # 500 bp sequence below --min-len
    assert t["plasme"]["precision_bp"] == 1.0
    assert t["union"]["recall_bp_assembly"] == 1.0


def test_unassigned_counts_as_false_positive(tmp_path):
    build(tmp_path)
    (tmp_path / "pred.paf").write_text(paf(("genomad|chim", 4000, 0, 1000, "PL1", 4000, 0, 1000)))
    evaluate.main([
        "--refs", str(tmp_path / "refs.fna"), "--provenance", str(tmp_path / "prov.tsv"),
        "--contigs-paf", str(tmp_path / "contigs.paf"), "--pred-paf", str(tmp_path / "pred.paf"),
        "--out-tsv", str(tmp_path / "m.tsv"),
    ])
    rows = (tmp_path / "m.tsv").read_text().splitlines()
    genomad = dict(zip(rows[0].split("\t"), rows[1].split("\t")))
    assert genomad["n_plasmid"] == "0"
    assert genomad["precision_bp"] == "0.0"


def test_merged_length():
    assert evaluate.merged_length([(0, 10), (5, 20), (30, 40)]) == 30
    assert evaluate.merged_length([]) == 0


def test_consensus_rows(tmp_path):
    build(tmp_path)
    # genomad|c1 and plasme|c2 do not overlap; add a third tool agreeing with genomad|c1 only
    (tmp_path / "pred.paf").write_text(paf(
        ("genomad|c1", 2500, 0, 2500, "PL1", 4000, 0, 2500),
        ("mobrecon|r1", 2400, 0, 2400, "PL1", 4000, 50, 2450),
        ("genomad|c3", 5000, 0, 5000, "CHR", 10000, 0, 5000),
        ("plasme|c2", 1500, 0, 1500, "PL1", 4000, 2500, 4000),
    ))
    # all-vs-all: genomad|c1 ~ mobrecon|r1 agree (>= 50 % of the shorter); a short hit c3 ~ c2 does not
    (tmp_path / "self.paf").write_text(paf(
        ("genomad|c1", 2500, 0, 2400, "mobrecon|r1", 2400, 0, 2400),
        ("genomad|c3", 5000, 0, 300, "plasme|c2", 1500, 0, 300),
    ))
    out = tmp_path / "m.json"
    evaluate.main([
        "--refs", str(tmp_path / "refs.fna"), "--provenance", str(tmp_path / "prov.tsv"),
        "--contigs-paf", str(tmp_path / "contigs.paf"), "--pred-paf", str(tmp_path / "pred.paf"),
        "--self-paf", str(tmp_path / "self.paf"),
        "--out-tsv", str(tmp_path / "m.tsv"), "--out-json", str(out),
    ])
    t = {r["tool"]: r for r in json.loads(out.read_text())["tools"]}
    assert t["min2"]["n_predicted"] == 2          # genomad|c1 + mobrecon|r1
    assert t["min2"]["precision_bp"] == 1.0       # the chromosome false positive is gone
    assert t["min3"]["n_predicted"] == 0
    assert t["union"]["n_predicted"] == 4


def test_provenance_with_host_column(tmp_path):
    (tmp_path / "prov.tsv").write_text(
        "type\taccession\tset\thost\tsource\nhost\tCHR\tt\t\tCHR\nplasmid\tPL1\tt\tCHR\tPL1\nplasmid\tPL2\tt\tCHR\tPL2\n")
    types = evaluate.read_provenance(str(tmp_path / "prov.tsv"))
    assert types == {"CHR": "host", "PL1": "plasmid", "PL2": "plasmid"}
