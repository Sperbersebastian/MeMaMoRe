#!/usr/bin/env python3
"""Smoke test for the Flask GUI against an empty, throwaway ROOT.

Run: python3 tests/gui_smoke.py   (needs flask)
"""
import os
import sys
import tempfile
from pathlib import Path

tmp = tempfile.mkdtemp(prefix="memamo_ci_")
os.environ["MEMAMO_ROOT"] = tmp
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "gui"))

import app as gui  # noqa: E402

c = gui.app.test_client()


def check(method, url, want, **kw):
    r = c.open(url, method=method, **kw)
    assert r.status_code == want, f"{method} {url}: {r.status_code} != {want}"


for page in ("/", "/samples", "/run", "/results", "/logs"):
    check("GET", page, 200)
check("GET", "/api/samples", 200)
check("GET", "/api/modules", 200)
check("GET", "/api/jobs", 200)

# Sample-name validation must reject traversal and odd characters
os.makedirs(os.path.join(tmp, "SRA", "qc", "fastp", "keep"))
for bad in ("..", ".hidden", "a b", "x;rm"):
    check("DELETE", f"/api/samples/{bad}", 400)
    check("GET", f"/api/status/{bad}", 400)
assert os.path.isdir(os.path.join(tmp, "SRA", "qc", "fastp", "keep"))

# Unknown module / bad input to /api/run is refused without starting a job
check("POST", "/api/run", 400, json={"sample": "mini", "module": "rm -rf"})
check("POST", "/api/run", 400, json={"sample": "../x", "module": "qc"})
assert not gui.JOBS, "a job was started for invalid input"

print("gui smoke: ok")
