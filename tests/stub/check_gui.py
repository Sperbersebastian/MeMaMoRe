#!/usr/bin/env python3
"""After tests/stub/run.sh: check the GUI API shows results for every module.

Usage: python3 tests/stub/check_gui.py <STUB_WORK> [sample]   (needs flask)
"""
import os
import sys
from pathlib import Path

work = sys.argv[1]
sample = sys.argv[2] if len(sys.argv) > 2 else "stub"
os.environ["MEMAMO_ROOT"] = work
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "gui"))
import app as gui  # noqa: E402

c = gui.app.test_client()
fail = []

samples = c.get("/api/samples").get_json()
names = [s["name"] for s in samples]
if sample not in names:
    fail.append(f"/api/samples does not list {sample!r}: {names}")
else:
    status = next(s["status"] for s in samples if s["name"] == sample)
    pending = [m for m, st in status.items() if st != "done"]
    if pending:
        fail.append(f"modules not 'done' in /api/samples: {pending} ({status})")

# result keys that must be non-empty for each module
expect = {
    "qc": ["fastp"],
    "assembly": ["assembly", "contigs_kept", "mapping"],
    "binning": ["magscot", "checkm2"],
    "plasmids": ["mobrecon", "genomad"],
    "viruses": ["checkv", "genomad"],
    "args": ["summary"],
}
for mod, keys in expect.items():
    r = c.get(f"/api/results/{sample}/{mod}")
    data = r.get_json() or {}
    for k in keys:
        if r.status_code != 200 or not data.get(k):
            fail.append(f"/api/results/{sample}/{mod}: '{k}' empty")
    print(f"{mod:9s} " + ", ".join(f"{k}={len(data.get(k) or [])}" if isinstance(data.get(k), list)
                                    else f"{k}={'ok' if data.get(k) else '-'}" for k in keys))

if fail:
    print("\n".join("FAIL " + f for f in fail))
    sys.exit(1)
print("stub gui check: ok")
