"""Unit tests for scripts/sim_plan.py (simulation config -> tables)."""
import csv
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("sim_plan", ROOT / "scripts" / "sim_plan.py")
sim_plan = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sim_plan)


def plan(tmp_path, text):
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(text)
    out = tmp_path / "plan"
    sim_plan.main(["sim_plan.py", str(cfg), str(out)])
    read = lambda n: list(csv.DictReader(open(out / n), delimiter="\t"))  # noqa: E731
    return read("references.tsv"), read("weights.tsv"), out


def test_repo_configs_parse(tmp_path):
    for cfg in [ROOT / "config" / "sim.yaml", *sorted((ROOT / "config" / "testsets").glob("*.yaml"))]:
        refs, weights, _ = plan(tmp_path, cfg.read_text())
        assert refs and weights, cfg


def test_legacy_single_sample(tmp_path):
    refs, weights, out = plan(tmp_path, (ROOT / "config" / "sim.yaml").read_text())
    assert {w["sample"] for w in weights} == {"miniset"}
    assert [r["type"] for r in refs].count("host") == 2
    assert all(r["host"] == "" for r in refs)
    assert "seed\t42" in (out / "settings.tsv").read_text()


def test_empty_hosts_mapping(tmp_path):
    refs, _, _ = plan(tmp_path, "set_name: x\nread_len: 150\nn_reads: 10\nhosts: {}\nplasmids: [P1]\n"
                                "abundances: {P1: 1.0}\n")
    assert refs == [{"accession": "P1", "type": "plasmid", "host": ""}]


EXT = """set_name: g
read_len: 150
n_reads: 1000
seed: 3
references:
  - {acc: H1, type: host, abundance: 0.6}
  - {acc: H2, type: host, abundance: 0.3}
  - {acc: P1, type: plasmid, host: H1, abundance: 0.1}
profiles: {n: 4, sigma: 1.0, copy_sigma: 0.0}
"""


def test_profiles_are_seeded_and_coupled(tmp_path):
    _, w1, _ = plan(tmp_path, EXT)
    _, w2, _ = plan(tmp_path, EXT)
    assert w1 == w2
    by = {}
    for r in w1:
        by.setdefault(r["sample"], {})[r["accession"]] = float(r["weight"])
    assert sorted(by) == ["g_s1", "g_s2", "g_s3", "g_s4"]
    for s in by.values():
        assert sum(s.values()) == pytest.approx(1.0, abs=1e-5)
        assert s["P1"] / s["H1"] == pytest.approx(0.1 / 0.6, rel=1e-3)  # copy_sigma 0: fixed ratio
    assert len({round(s["H1"], 4) for s in by.values()}) > 1


def test_unknown_host_is_rejected(tmp_path):
    with pytest.raises(SystemExit):
        plan(tmp_path, EXT.replace("host: H1", "host: H9"))
