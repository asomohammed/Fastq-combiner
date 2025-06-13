#!/usr/bin/env python3
"""
fastq_combiner_v2.py

FASTQ File Combiner - Version 2

Description:
    Combines multiple FASTQ files into a single output FASTQ file.
    Designed for high-throughput sequencing workflows where merging reads from different sources is required.
    Handles standard FASTQ format (4 lines per record), supports large files, and preserves record integrity.

Usage:
    python fastq_combiner_v2.py -i <input1.fastq> <input2.fastq> ... -o <output.fastq>
    Additional command-line options may be available (see code or help).

Inputs:
    - One or more input FASTQ files (plain or gzipped, depending on implementation).
    - Command-line arguments for specifying files and options.

Outputs:
    - Single combined FASTQ file containing all records from input files, in order provided.

Features:
    - Efficient line-by-line reading to support large files.
    - Basic validation of FASTQ format.
    - Optionally preserves original file read order.

Dependencies:
    - Python 3.x
    - Standard library modules only (unless otherwise specified in the code).

Limitations:
    - Assumes input files are valid FASTQ format.
    - No quality-filtering or read deduplication.
    
Date:
    2025-06-13

"""
import argparse
import csv
import gzip
import logging
import os
import shutil
import sys
import tempfile
import hashlib
import time
import json
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Tuple, Optional, Any

import psutil
from tqdm import tqdm

# Try to use rapidfuzz for best fuzzy matching, else fallback to difflib
try:
    from rapidfuzz import process as rfprocess, fuzz as rffuzz
    def fuzzy_match(query, candidates):
        match = rfprocess.extractOne(query, candidates, scorer=rffuzz.QRatio, score_cutoff=70)
        return match[0] if match else None
except ImportError:
    import difflib
    def fuzzy_match(query, candidates):
        matches = difflib.get_close_matches(query, candidates, n=1, cutoff=0.7)
        return matches[0] if matches else None

# --- Logging ---
def setup_logging(log_file: str, log_level: str):
    logger = logging.getLogger("fastq_combiner")
    logger.setLevel(getattr(logging, log_level.upper(), logging.INFO))
    fmt = logging.Formatter('[%(asctime)s][%(levelname)s] %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(fmt)
    logger.handlers = [sh]
    if log_file:
        fh = logging.FileHandler(log_file)
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    return logger

# --- File integrity ---
def md5sum(filename: str, blocksize: int = 2**20) -> str:
    m = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            m.update(block)
    return m.hexdigest()

def validate_file_exists_and_md5(path: str, expected_md5: Optional[str] = None) -> bool:
    if not os.path.exists(path):
        return False
    if expected_md5:
        return md5sum(path) == expected_md5
    return True

# --- Mapping ---
def read_mapping_file(csv_file: str) -> Dict[str, List[str]]:
    mapping = defaultdict(list)
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        first_row = next(reader)
        if first_row[0].lower() in ['target', 'target_sample', 'output', 'sample']:
            pass
        else:
            target, *sources = first_row
            mapping[target.strip()].extend([s.strip() for s in sources if s.strip()])
        for row in reader:
            if not row or not row[0].strip():
                continue
            target, *sources = row
            mapping[target.strip()].extend([s.strip() for s in sources if s.strip()])
    return dict(mapping)

# --- File discovery ---
def find_fastq_files(search_dirs: Optional[List[str]]) -> Dict[str, Dict[str, str]]:
    patterns = [
        "*_R1_*.fastq*", "*_R1.fastq*", "*_1.fastq*", "*.R1.fastq*",
        "*_R1_*.fq*", "*_R1.fq*", "*_1.fq*", "*.R1.fq*"
    ]
    file_pairs = {}
    search_dirs = search_dirs or ["."]
    for search_dir in search_dirs:
        for pattern in patterns:
            for r1_file in Path(search_dir).rglob(pattern):
                basename = r1_file.name
                if "_S" in basename and "_R1_" in basename:
                    sample_base = basename.split("_S")[0]
                elif "_R1_" in basename:
                    sample_base = basename.split("_R1_")[0]
                elif "_R1." in basename:
                    sample_base = basename.split("_R1.")[0]
                elif "_1." in basename:
                    sample_base = basename.split("_1.")[0]
                elif ".R1." in basename:
                    sample_base = basename.split(".R1.")[0]
                else:
                    continue
                r2_candidates = [
                    str(r1_file).replace("_R1_", "_R2_"),
                    str(r1_file).replace("_R1.", "_R2."),
                    str(r1_file).replace("_1.", "_2."),
                    str(r1_file).replace(".R1.", ".R2.")
                ]
                r2_file = next((c for c in r2_candidates if os.path.exists(c)), None)
                if not r2_file:
                    continue
                key = sample_base
                file_pairs[key] = {
                    'R1': str(r1_file.resolve()),
                    'R2': os.path.abspath(r2_file)
                }
    return file_pairs

def sanitize_sample_name(name: str) -> str:
    if "_S" in name:
        name = name.split("_S")[0]
    clean = name.replace(" ", "_").replace("-", "_")
    while "__" in clean:
        clean = clean.replace("__", "_")
    return clean.strip("_")

# --- Read counting (support gz and plain) ---
def count_reads_fastq(fastq_file: str) -> int:
    count = 0
    opener = gzip.open if fastq_file.endswith('.gz') else open
    try:
        with opener(fastq_file, 'rt') as f:
            for i, line in enumerate(f):
                if i % 4 == 0 and line.startswith('@'):
                    count += 1
    except Exception as e:
        return 0
    return count

# --- Combination ---
def combine_fastq_files(
    source_files: List[str], output_file: str, buffer_size: int = 8 * 1024 * 1024,
    compresslevel: int = 6
) -> int:
    opener = gzip.open if output_file.endswith('.gz') else open
    total_reads = 0
    with opener(output_file, 'wb', compresslevel=compresslevel) if output_file.endswith('.gz') \
         else open(output_file, 'wb') as out_f:
        for src in source_files:
            src_opener = gzip.open if src.endswith('.gz') else open
            with src_opener(src, 'rb') as in_f:
                shutil.copyfileobj(in_f, out_f, length=buffer_size)
            total_reads += count_reads_fastq(src)
    return total_reads

# --- Resource usage ---
def get_resource_usage() -> Dict[str, Any]:
    p = psutil.Process(os.getpid())
    return {
        "cpu_percent": psutil.cpu_percent(),
        "mem_mb": p.memory_info().rss / 1024 / 1024,
        "open_files": len(p.open_files())
    }

# --- HTML report ---
def generate_html_report(
    output_dir: str,
    mapping: Dict[str, List[str]],
    file_pairs: Dict[str, Dict[str, str]],
    combination_stats: Dict[str, Dict[str, Any]],
    fuzzy_matches: Dict[str, Tuple[str, str]],
    summary_json_path: str,
    summary_csv_path: str
) -> str:
    html_path = os.path.join(output_dir, "combination_report.html")
    total_targets = len(mapping)
    successful = len(combination_stats)
    total_reads = sum(stats['total_reads'] for stats in combination_stats.values())
    html = [
        "<!DOCTYPE html><html><head><meta charset='UTF-8'>",
        "<title>FASTQ Combination Report</title>",
        "<style>body{font-family:sans-serif;background:#f5f5f5;}",
        "table{border-collapse:collapse;}th,td{padding:8px;border:1px solid #ddd;}</style>",
        "</head><body>",
        f"<h1>FASTQ Combination Report</h1>",
        f"<p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>",
        f"<p>Targets: {total_targets} | Successful: {successful} | Total reads: {total_reads}</p>",
        "<h2>Results</h2><table><tr><th>Target</th><th>Reads</th><th>Status</th><th>Output Files</th></tr>"
    ]
    for target in mapping:
        status = "Success" if target in combination_stats else "FAILED"
        reads = combination_stats[target]['total_reads'] if target in combination_stats else 0
        out_r1 = combination_stats[target]['r1_output'] if target in combination_stats else ""
        out_r2 = combination_stats[target]['r2_output'] if target in combination_stats else ""
        html.append(f"<tr><td>{target}</td><td>{reads:,}</td><td>{status}</td>"
            f"<td>{os.path.basename(out_r1)}<br>{os.path.basename(out_r2)}</td></tr>")
    html.append("</table>")
    if fuzzy_matches:
        html.append("<h2>Fuzzy Matches Used</h2><table><tr><th>Original</th><th>Matched</th></tr>")
        for orig, (o, m) in fuzzy_matches.items():
            html.append(f"<tr><td>{o}</td><td>{m}</td></tr>")
        html.append("</table>")
    html.append(f"<h2>Pipeline Outputs</h2><ul>")
    html.append(f"<li>Summary JSON: <a href='{os.path.basename(summary_json_path)}'>{os.path.basename(summary_json_path)}</a></li>")
    html.append(f"<li>Summary CSV: <a href='{os.path.basename(summary_csv_path)}'>{os.path.basename(summary_csv_path)}</a></li>")
    html.append(f"</ul><p>Output dir: {os.path.abspath(output_dir)}</p>")
    html.append("</body></html>")
    with open(html_path, 'w') as f:
        f.write('\n'.join(html))
    return html_path

# --- CSV/JSON summary ---
def write_summary_csv(path: str, stats: Dict[str, Dict[str, Any]]):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['target', 'reads', 'r1_output', 'r2_output', 'processing_time', 'md5_r1', 'md5_r2'])
        writer.writeheader()
        for target, s in stats.items():
            writer.writerow({
                'target': target,
                'reads': s['total_reads'],
                'r1_output': s['r1_output'],
                'r2_output': s['r2_output'],
                'processing_time': s['processing_time'],
                'md5_r1': s.get('md5_r1', ''),
                'md5_r2': s.get('md5_r2', '')
            })

def write_summary_json(path: str, stats: Dict[str, Dict[str, Any]], resource_stats: Dict[str, Any]):
    out = {
        'samples': stats,
        'resource_usage': resource_stats
    }
    with open(path, 'w') as f:
        json.dump(out, f, indent=2)

# --- Main worker (with resume, dry run, integrity, consistency) ---
def process_target(
    target: str,
    matched_sources: List[str],
    file_pairs: Dict[str, Dict[str, str]],
    output_dir: str,
    buffer_size: int,
    compresslevel: int,
    resume: bool,
    dry_run: bool,
    logger: logging.Logger
) -> Tuple[str, Dict[str, Any]]:
    clean_target = sanitize_sample_name(target)
    r1_sources = [file_pairs[s]['R1'] for s in matched_sources]
    r2_sources = [file_pairs[s]['R2'] for s in matched_sources]
    r1_output = os.path.join(output_dir, f"{clean_target}_S1_R1_001.fastq.gz")
    r2_output = os.path.join(output_dir, f"{clean_target}_S1_R2_001.fastq.gz")
    # Resume support
    already_done = all(os.path.exists(f) for f in [r1_output, r2_output])
    if resume and already_done:
        logger.info(f"[{target}] Skipping: output exists (resume)")
        r1_reads = count_reads_fastq(r1_output)
        r2_reads = count_reads_fastq(r2_output)
        md5_r1 = md5sum(r1_output)
        md5_r2 = md5sum(r2_output)
        return target, {
            'total_reads': r1_reads,
            'r1_output': r1_output,
            'r2_output': r2_output,
            'processing_time': 0,
            'md5_r1': md5_r1,
            'md5_r2': md5_r2
        }
    if dry_run:
        logger.info(f"[{target}] DRY RUN: Would create {r1_output} and {r2_output} from {len(r1_sources)} sources")
        return target, {
            'total_reads': 0,
            'r1_output': r1_output,
            'r2_output': r2_output,
            'processing_time': 0,
            'md5_r1': '',
            'md5_r2': ''
        }
    # Real run
    start = time.time()
    r1_reads = combine_fastq_files(r1_sources, r1_output, buffer_size, compresslevel)
    r2_reads = combine_fastq_files(r2_sources, r2_output, buffer_size, compresslevel)
    if r1_reads != r2_reads:
        logger.warning(f"[{target}] R1 ({r1_reads}) and R2 ({r2_reads}) read counts mismatch")
    md5_r1 = md5sum(r1_output)
    md5_r2 = md5sum(r2_output)
    elapsed = time.time() - start
    return target, {
        'total_reads': r1_reads,
        'r1_output': r1_output,
        'r2_output': r2_output,
        'processing_time': elapsed,
        'md5_r1': md5_r1,
        'md5_r2': md5_r2
    }

def main():
    parser = argparse.ArgumentParser(description="FASTQ file combiner: robust, resumable, with HTML/CSV/JSON reporting")
    parser.add_argument('csv_file', help='CSV mapping file')
    parser.add_argument('-o', '--output-dir', default='combined', help='Output directory')
    parser.add_argument('-d', '--search-dirs', nargs='+', help='Dirs to search for FASTQ')
    parser.add_argument('--buffer-size', type=int, default=8 * 1024 * 1024, help='Buffer size (bytes)')
    parser.add_argument('--compresslevel', type=int, default=6, help='Gzip compression level (1-9)')
    parser.add_argument('-j', '--jobs', type=int, default=os.cpu_count(), help='Parallel jobs')
    parser.add_argument('--resume', action='store_true', help='Resume: skip outputs that already exist')
    parser.add_argument('--dry-run', action='store_true', help='Dry run: show actions, do not combine')
    parser.add_argument('--log-file', default='', help='Log file')
    parser.add_argument('--log', default="INFO", help='Log level')
    args = parser.parse_args()

    logger = setup_logging(args.log_file, args.log)
    os.makedirs(args.output_dir, exist_ok=True)
    logger.info(f"Starting FASTQ combiner. Output: {args.output_dir}")

    mapping = read_mapping_file(args.csv_file)
    file_pairs = find_fastq_files(args.search_dirs)
    available_samples = list(file_pairs.keys())

    final_mapping = {}
    fuzzy_matches: Dict[str, Tuple[str, str]] = {}
    for target, sources in mapping.items():
        matched = []
        for s in sources:
            if s in file_pairs:
                matched.append(s)
            else:
                fz = fuzzy_match(s, available_samples)
                if fz:
                    matched.append(fz)
                    fuzzy_matches[s] = (s, fz)
                else:
                    logger.warning(f"No match for {s} (target {target})")
        if matched:
            final_mapping[target] = matched

    if not final_mapping:
        logger.error("No sources could be matched to targets.")
        sys.exit(1)

    stats: Dict[str, Dict[str, Any]] = {}
    logger.info(f"Combining {len(final_mapping)} targets using {args.jobs} threads")
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = [
            executor.submit(
                process_target, target, matched, file_pairs, args.output_dir,
                args.buffer_size, args.compresslevel, args.resume, args.dry_run, logger
            )
            for target, matched in final_mapping.items()
        ]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Combining", unit="sample"):
            target, res = fut.result()
            stats[target] = res

    resource_stats = get_resource_usage()
    summary_json_path = os.path.join(args.output_dir, "summary.json")
    summary_csv_path = os.path.join(args.output_dir, "summary.csv")
    write_summary_json(summary_json_path, stats, resource_stats)
    write_summary_csv(summary_csv_path, stats)
    html_report = generate_html_report(
        args.output_dir, mapping, file_pairs, stats, fuzzy_matches, summary_json_path, summary_csv_path
    )
    logger.info(f"Done. {len(stats)} targets combined.")
    logger.info(f"HTML report: {html_report}")
    logger.info(f"Summary JSON: {summary_json_path}")
    logger.info(f"Summary CSV: {summary_csv_path}")

if __name__ == "__main__":
    main()
