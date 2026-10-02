#!/usr/bin/env python3
"""MeMaMoRe Web GUI – Flask backend."""
import json, os, re, subprocess, threading, time, glob, csv, signal, shutil
from pathlib import Path
from flask import Flask, render_template, jsonify, request, Response, abort

app = Flask(__name__)
ROOT = os.environ.get("MEMAMO_ROOT", str(Path(__file__).resolve().parent.parent))

# ── running jobs ──────────────────────────────────────────────────────────────
JOBS = {}  # id -> {proc, sample, module, start, log}
_job_id = 0
_lock = threading.Lock()

MODULES = [
    {"id": "ingest",   "name": "Ingest",   "desc": "Import & validate reads"},
    {"id": "qc",       "name": "QC",       "desc": "Quality control & trimming"},
    {"id": "assembly", "name": "Assembly", "desc": "metaSPAdes assembly"},
    {"id": "binning",  "name": "Binning",  "desc": "MetaBat2 + MAGScoT"},
    {"id": "plasmids", "name": "Plasmids", "desc": "Plasmid detection"},
    {"id": "viruses",  "name": "Viruses",  "desc": "Virus detection"},
    {"id": "args",     "name": "ARGs",     "desc": "Resistance gene detection"},
]

MODULE_IDS = {m["id"] for m in MODULES}

# ── input validation ──────────────────────────────────────────────────────────
# Sample names end up in filesystem paths (and rmtree), so only allow a safe
# charset; the leading alnum rules out "." and "..".
_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$")

def _valid_name(name):
    return isinstance(name, str) and bool(_NAME_RE.match(name))

def _safe_sra_path(*parts):
    """Join parts under ROOT/SRA and refuse anything that escapes it."""
    base = os.path.realpath(os.path.join(ROOT, "SRA"))
    path = os.path.realpath(os.path.join(base, *parts))
    if os.path.commonpath([base, path]) != base or path == base:
        abort(400, description="Invalid path")
    return path

@app.before_request
def _validate_path_args():
    for key in ("sample", "name"):
        val = (request.view_args or {}).get(key)
        if val is not None and not _valid_name(val):
            abort(400, description=f"Invalid {key}")

@app.errorhandler(400)
def _bad_request(e):
    return jsonify({"error": e.description}), 400

# ── helpers ───────────────────────────────────────────────────────────────────
def _samples():
    """Discover samples from SRA subdirectories."""
    samples = set()
    for sub in ["assemblies/spades", "qc/fastp", "viruses", "args", "plasmids"]:
        d = os.path.join(ROOT, "SRA", sub)
        if os.path.isdir(d):
            for name in os.listdir(d):
                if os.path.isdir(os.path.join(d, name)) and not name.startswith("."):
                    samples.add(name)
    return sorted(samples)

def _module_status(sample):
    """Check .done sentinels for each module's substeps."""
    status = {}
    for mod in MODULES:
        mid = mod["id"]
        # Check various done-file locations
        done_paths = [
            os.path.join(ROOT, "SRA", mid, sample, ".done"),
            os.path.join(ROOT, "SRA", mid, sample),
        ]
        # Module-specific done files
        if mid == "qc":
            done = os.path.join(ROOT, "SRA", "qc", "fastp", sample)
            status[mid] = "done" if os.path.isdir(done) and os.listdir(done) else "pending"
        elif mid == "assembly":
            contigs = os.path.join(ROOT, "SRA", "assemblies", "spades", sample, "contigs.fasta")
            status[mid] = "done" if os.path.isfile(contigs) else "pending"
        elif mid == "binning":
            bins = glob.glob(os.path.join(ROOT, "SRA", "binning", "metabat2", sample, "bin.*.fa"))
            status[mid] = "done" if bins else "pending"
        elif mid in ("viruses", "args", "plasmids"):
            done_dir = os.path.join(ROOT, "SRA", mid, sample, ".done")
            if os.path.isdir(done_dir) and os.listdir(done_dir):
                status[mid] = "done"
            else:
                status[mid] = "pending"
        elif mid == "ingest":
            # Check if any reads exist
            qc_dir = os.path.join(ROOT, "SRA", "qc", "fastp", sample)
            raw_dir = os.path.join(ROOT, "SRA", "reads", sample)
            if os.path.isdir(qc_dir) or os.path.isdir(raw_dir):
                status[mid] = "done"
            else:
                status[mid] = "pending"
        else:
            status[mid] = "pending"

        # Check if currently running
        with _lock:
            for j in JOBS.values():
                if j["sample"] == sample and j["module"] == mid and j["proc"].poll() is None:
                    status[mid] = "running"
    return status

def _read_tsv(path, max_rows=5000):
    """Read a TSV file and return list of dicts."""
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for i, row in enumerate(reader):
            if i >= max_rows:
                break
            rows.append(dict(row))
    return rows

# ── pages ─────────────────────────────────────────────────────────────────────
@app.route("/")
def dashboard():
    return render_template("dashboard.html")

@app.route("/samples")
def samples_page():
    return render_template("samples.html")

@app.route("/run")
def run_page():
    return render_template("run.html")

@app.route("/results")
def results_page():
    return render_template("results.html")

@app.route("/logs")
def logs_page():
    return render_template("logs.html")

# ── API ───────────────────────────────────────────────────────────────────────
@app.route("/api/samples")
def api_samples():
    samples = _samples()
    result = []
    for s in samples:
        result.append({"name": s, "status": _module_status(s)})
    return jsonify(result)

@app.route("/api/modules")
def api_modules():
    return jsonify(MODULES)

@app.route("/api/status/<sample>")
def api_status(sample):
    return jsonify(_module_status(sample))

@app.route("/api/run", methods=["POST"])
def api_run():
    global _job_id
    data = request.get_json(silent=True) or {}
    sample = data.get("sample", "")
    samples_list = data.get("samples", [])  # multi-sample support
    modules = data.get("modules", [])
    force = data.get("force", False)
    cpus = data.get("cpus", 14)

    if not sample and not samples_list:
        return jsonify({"error": "No sample specified"}), 400
    if not modules:
        return jsonify({"error": "No modules selected"}), 400
    if not isinstance(modules, list) or any(m not in MODULE_IDS for m in modules):
        return jsonify({"error": "Unknown module"}), 400
    if sample and sample != "ALL" and not _valid_name(sample):
        return jsonify({"error": "Invalid sample name"}), 400
    if not isinstance(samples_list, list) or not all(_valid_name(s) for s in samples_list):
        return jsonify({"error": "Invalid sample name"}), 400
    try:
        cpus = int(cpus)
    except (TypeError, ValueError):
        return jsonify({"error": "cpus must be an integer"}), 400
    if not 1 <= cpus <= 1024:
        return jsonify({"error": "cpus out of range"}), 400

    if samples_list:
        target_samples = samples_list
    else:
        target_samples = _samples() if sample == "ALL" else [sample]

    # Check dependencies and collect warnings
    warnings = []
    for s in target_samples:
        status = _module_status(s)
        for mod in modules:
            deps = MODULE_DEPS.get(mod, [])
            for dep in deps:
                if status.get(dep) != "done":
                    warnings.append(f"{s}: '{mod}' requires '{dep}' which is not done yet")

    launched = []

    for s in target_samples:
        for mod in modules:
            cmd = [os.path.join(ROOT, "bin", "main.sh"), "run", mod, "--sample", s, "--cpus", str(cpus)]
            if force:
                cmd.append("--force")

            log_dir = os.path.join(ROOT, "logs")
            os.makedirs(log_dir, exist_ok=True)
            log_path = os.path.join(log_dir, f"gui_{mod}_{s}.log")

            with open(log_path, "w") as lf:
                proc = subprocess.Popen(
                    cmd, stdout=lf, stderr=subprocess.STDOUT,
                    cwd=ROOT, preexec_fn=os.setsid
                )

            with _lock:
                _job_id += 1
                jid = _job_id
                JOBS[jid] = {
                    "proc": proc, "sample": s, "module": mod,
                    "start": time.time(), "log": log_path, "pid": proc.pid
                }
            launched.append({"id": jid, "module": mod, "sample": s})

    result = {"launched": launched}
    if warnings:
        result["warnings"] = warnings
    return jsonify(result)

@app.route("/api/jobs")
def api_jobs():
    result = []
    with _lock:
        for jid, j in JOBS.items():
            poll = j["proc"].poll()
            result.append({
                "id": jid,
                "sample": j["sample"],
                "module": j["module"],
                "running": poll is None,
                "exit_code": poll,
                "elapsed": int(time.time() - j["start"]),
                "log": j["log"],
            })
    return jsonify(result)

@app.route("/api/jobs/<int:jid>", methods=["DELETE"])
def api_cancel_job(jid):
    with _lock:
        j = JOBS.get(jid)
    if not j:
        return jsonify({"error": "Job not found"}), 404
    try:
        os.killpg(os.getpgid(j["proc"].pid), signal.SIGTERM)
    except ProcessLookupError:
        pass
    return jsonify({"cancelled": jid})

@app.route("/api/results/<sample>/<module>")
def api_results(sample, module):
    """Return results for a given sample and module."""
    results = {}
    if module == "qc":
        # Read fastp JSON stats
        fastp_json = os.path.join(ROOT, "SRA", "qc", "fastp", sample, f"{sample}.fastp.json")
        if not os.path.isfile(fastp_json):
            # Try alternate naming
            candidates = glob.glob(os.path.join(ROOT, "SRA", "qc", "fastp", sample, "*.fastp.json"))
            if candidates:
                fastp_json = candidates[0]
        if os.path.isfile(fastp_json):
            try:
                with open(fastp_json) as f:
                    fj = json.load(f)
                s = fj.get("summary", {})
                bf, af = s.get("before_filtering", {}), s.get("after_filtering", {})
                results["fastp"] = {
                    "total_reads_before": bf.get("total_reads", 0),
                    "total_reads_after": af.get("total_reads", 0),
                    "total_bases_before": bf.get("total_bases", 0),
                    "total_bases_after": af.get("total_bases", 0),
                    "q30_before": bf.get("q30_rate", 0),
                    "q30_after": af.get("q30_rate", 0),
                    "gc_before": bf.get("gc_content", 0),
                    "gc_after": af.get("gc_content", 0),
                    "duplication_rate": fj.get("duplication", {}).get("rate", 0),
                    "adapter_trimmed_reads": fj.get("adapter_cutting", {}).get("adapter_trimmed_reads", 0),
                    "reads_passed_filter": fj.get("filtering_result", {}).get("passed_filter_reads", 0),
                    "reads_low_quality": fj.get("filtering_result", {}).get("low_quality_reads", 0),
                    "reads_too_short": fj.get("filtering_result", {}).get("too_short_reads", 0),
                }
            except (json.JSONDecodeError, KeyError):
                results["fastp"] = None
        else:
            results["fastp"] = None
    elif module == "assembly":
        contigs_path = os.path.join(ROOT, "SRA", "assemblies", "spades", sample, "contigs.fasta")
        if os.path.isfile(contigs_path):
            lengths = []
            with open(contigs_path) as f:
                cur_len = 0
                for line in f:
                    if line.startswith(">"):
                        if cur_len > 0:
                            lengths.append(cur_len)
                        cur_len = 0
                    else:
                        cur_len += len(line.strip())
                if cur_len > 0:
                    lengths.append(cur_len)
            lengths.sort(reverse=True)
            total = sum(lengths)
            # N50
            n50 = 0
            cumsum = 0
            for l in lengths:
                cumsum += l
                if cumsum >= total / 2:
                    n50 = l
                    break
            results["assembly"] = {
                "num_contigs": len(lengths),
                "total_length": total,
                "n50": n50,
                "longest": lengths[0] if lengths else 0,
                "shortest": lengths[-1] if lengths else 0,
            }
        else:
            results["assembly"] = None
        # Filtered contigs
        filtered_path = os.path.join(ROOT, "SRA", "assemblies", "contig_qc", sample, "contigs_kept.tsv")
        results["contigs_kept"] = _read_tsv(filtered_path)
        mapping_path = os.path.join(ROOT, "SRA", "assemblies", "contig_qc", sample, "mapping_counts.tsv")
        results["mapping"] = _read_tsv(mapping_path)
    elif module == "args":
        summary_path = os.path.join(ROOT, "SRA", "args", sample, "summary", "args_summary.tsv")
        results["summary"] = _read_tsv(summary_path)
        # Per tool
        for tool in ["deeparg", "rgi", "amrplusplus"]:
            tool_dir = os.path.join(ROOT, "SRA", "args", sample, tool)
            if os.path.isdir(tool_dir):
                for ctx in ["contigs", "mags", "plasmids", "viruses"]:
                    ctx_dir = os.path.join(tool_dir, ctx)
                    if os.path.isdir(ctx_dir):
                        for f in glob.glob(os.path.join(ctx_dir, "*.txt")) + glob.glob(os.path.join(ctx_dir, "*.mapping.ARG")):
                            key = f"{tool}/{ctx}/{Path(f).stem}"
                            results[key] = _read_tsv(f)
                # AMR++ specific
                for f in glob.glob(os.path.join(tool_dir, "*.tsv")):
                    key = f"{tool}/{Path(f).stem}"
                    results[key] = _read_tsv(f)
    elif module == "viruses":
        checkv_path = os.path.join(ROOT, "SRA", "viruses", sample, "checkv", "quality_summary.tsv")
        results["checkv"] = _read_tsv(checkv_path)
        genomad_path = os.path.join(ROOT, "SRA", "viruses", sample, "genomad", "contigs_summary", "contigs_virus_summary.tsv")
        results["genomad"] = _read_tsv(genomad_path)
    elif module == "plasmids":
        mobrecon_path = os.path.join(ROOT, "SRA", "plasmids", sample, "mobrecon", "mobtyper_results.txt")
        results["mobrecon"] = _read_tsv(mobrecon_path)
        # geNomad names outputs after the input file (contigs_summary/contigs_plasmid_summary.tsv)
        hits = glob.glob(os.path.join(ROOT, "SRA", "plasmids", sample, "genomad", "*_summary", "*_plasmid_summary.tsv"))
        results["genomad"] = _read_tsv(hits[0]) if hits else []
    elif module == "binning":
        magscot_path = os.path.join(ROOT, "SRA", "binning", "magscot", sample, "contigs_to_bin.tsv")
        results["magscot"] = _read_tsv(magscot_path)
        # CheckM2 quality stats
        checkm2_path = os.path.join(ROOT, "SRA", "binning", "qc_checkm2", sample, "quality_report.tsv")
        results["checkm2"] = _read_tsv(checkm2_path)
    return jsonify(results)

@app.route("/api/logs/<sample>")
def api_logs_stream(sample):
    """SSE stream of log file."""
    module = request.args.get("module", "")
    if module and module not in MODULE_IDS:
        return jsonify({"error": "Unknown module"}), 400
    if module:
        log_path = os.path.join(ROOT, "logs", f"gui_{module}_{sample}.log")
        if not os.path.isfile(log_path):
            log_path = os.path.join(ROOT, "logs", f"{module}_{sample}.log")
    else:
        log_path = os.path.join(ROOT, "logs", f"gui_*_{sample}.log")
        candidates = sorted(glob.glob(log_path), key=os.path.getmtime, reverse=True)
        log_path = candidates[0] if candidates else ""

    if not log_path or not os.path.isfile(log_path):
        return jsonify({"error": "No log file found"}), 404

    def generate():
        with open(log_path) as f:
            # Send existing content
            content = f.read()
            if content:
                yield f"data: {json.dumps({'lines': content.splitlines()[-200:]})}\n\n"
            # Then tail for new lines
            while True:
                line = f.readline()
                if line:
                    yield f"data: {json.dumps({'line': line.rstrip()})}\n\n"
                else:
                    time.sleep(0.5)
                    # Check if any job is still running for this sample
                    any_running = False
                    with _lock:
                        for j in JOBS.values():
                            if j["sample"] == sample and j["proc"].poll() is None:
                                any_running = True
                    if not any_running and not f.readline():
                        yield f"data: {json.dumps({'done': True})}\n\n"
                        break

    return Response(generate(), mimetype="text/event-stream",
                    headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"})

@app.route("/api/samples/<name>", methods=["DELETE"])
def api_delete_sample(name):
    """Delete a sample and all its output directories."""
    removed = []
    for sub in ["qc/fastp", "qc/fastqc", "assemblies/spades", "assemblies/contig_qc",
                "viruses", "args", "plasmids", "reads",
                "binning/comebin", "binning/binny", "binning/metabat2",
                "binning/magscot", "binning/qc_checkm2", "binning/coverage",
                "binning/gtdbtk", "binning/qc_rrna_trna"]:
        d = _safe_sra_path(sub, name)
        if os.path.isdir(d) and not os.path.islink(d):
            shutil.rmtree(d)
            removed.append(sub)
    return jsonify({"ok": True, "sample": name, "removed": removed})

@app.route("/api/disk/<sample>")
def api_disk(sample):
    """Return total disk usage in bytes for a sample."""
    total = 0
    for sub in ["qc/fastp", "qc/fastqc", "assemblies/spades", "assemblies/contig_qc",
                "viruses", "args", "plasmids", "reads",
                "binning/comebin", "binning/binny", "binning/metabat2",
                "binning/magscot", "binning/qc_checkm2"]:
        d = os.path.join(ROOT, "SRA", sub, sample)
        if os.path.isdir(d):
            for dirpath, dirnames, filenames in os.walk(d):
                for fn in filenames:
                    try:
                        total += os.path.getsize(os.path.join(dirpath, fn))
                    except OSError:
                        pass
    return jsonify({"sample": sample, "bytes": total})

# ── Pipeline dependency map ──────────────────────────────────────────────────
MODULE_DEPS = {
    "qc":       ["ingest"],
    "assembly": ["qc"],
    "binning":  ["assembly"],
    "plasmids": ["assembly"],
    "viruses":  ["assembly"],
    "args":     ["assembly"],
}

@app.route("/api/samples/add", methods=["POST"])
def api_add_sample():
    """Add a new sample via SRA accession, local path, or direct FASTA upload."""
    if request.content_type and "multipart/form-data" in request.content_type:
        sample_name = request.form.get("name", "").strip()
        source_type = request.form.get("type", "")
        file_obj = request.files.get("file")
        if not sample_name or not file_obj:
            return jsonify({"error": "Sample name and file required"}), 400
        if not _valid_name(sample_name):
            return jsonify({"error": "Invalid sample name (allowed: A-Z a-z 0-9 . _ -)"}), 400

        # Create base dirs
        for d in ["qc/fastp", "assemblies/spades", "viruses", "args", "plasmids"]:
            os.makedirs(os.path.join(ROOT, "SRA", d, sample_name), exist_ok=True)

        out_dir = os.path.join(ROOT, "SRA", "assemblies", "spades", sample_name)
        os.makedirs(out_dir, exist_ok=True) # Ensure output directory exists
        out_path = os.path.join(out_dir, "contigs.fasta")
        file_obj.save(out_path)

        return jsonify({"ok": True, "sample": sample_name, "note": "FASTA uploaded. Ready for Binning/Viruses."})

    data = request.get_json(silent=True) or {}
    sample_name = str(data.get("name", "")).strip()
    source_type = data.get("type", "")  # "sra" or "local"
    source_value = str(data.get("value", "")).strip()

    if not sample_name:
        return jsonify({"error": "Sample name required"}), 400
    if not _valid_name(sample_name):
        return jsonify({"error": "Invalid sample name (allowed: A-Z a-z 0-9 . _ -)"}), 400

    if source_type == "local":
        # Just create the sample directory structure
        for d in ["qc/fastp", "assemblies/spades", "viruses", "args", "plasmids"]:
            os.makedirs(os.path.join(ROOT, "SRA", d, sample_name), exist_ok=True)
        return jsonify({"ok": True, "sample": sample_name, "note": "Directories created. Run ingest to import data."})

    elif source_type == "sra":
        if not re.match(r"^[A-Za-z0-9]+$", source_value):
            return jsonify({"error": "Invalid SRA accession"}), 400
        # Launch ingest with --srr
        cmd = [os.path.join(ROOT, "bin", "main.sh"), "run", "ingest",
               "--sample", sample_name, "--srr", source_value]
        log_path = os.path.join(ROOT, "logs", f"gui_ingest_{sample_name}.log")
        os.makedirs(os.path.dirname(log_path), exist_ok=True)
        with open(log_path, "w") as lf:
            proc = subprocess.Popen(cmd, stdout=lf, stderr=subprocess.STDOUT,
                                    cwd=ROOT, preexec_fn=os.setsid)
        global _job_id
        with _lock:
            _job_id += 1
            jid = _job_id
            JOBS[jid] = {"proc": proc, "sample": sample_name, "module": "ingest",
                         "start": time.time(), "log": log_path, "pid": proc.pid}
        return jsonify({"ok": True, "sample": sample_name, "job_id": jid})

    return jsonify({"error": "Invalid source type"}), 400

# ── main ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--port", type=int, default=5000)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--root", default=ROOT)
    args = parser.parse_args()
    ROOT = args.root
    os.environ["MEMAMO_ROOT"] = ROOT
    print(f"[MeMaMoRe GUI] ROOT={ROOT}")
    if args.host not in ("127.0.0.1", "localhost", "::1"):
        print("[MeMaMoRe GUI] WARNING: listening on a non-local address. The GUI has no "
              "authentication; anyone who can reach this port can start jobs and delete samples.")
    print(f"[MeMaMoRe GUI] http://{args.host}:{args.port}")
    app.run(host=args.host, port=args.port, debug=False, threaded=True)
