import argparse
import os
import csv
from pathlib import Path
from .utils import count_reads_fastq, combine_fastq_files, md5sum
from .report import generate_html_report


def read_mapping_file(csv_file):
    mapping = {}
    with open(csv_file, newline='') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or not row[0].strip():
                continue
            target, *sources = row
            mapping[target.strip()] = [s.strip() for s in sources if s.strip()]
    return mapping

def find_fastq_pairs(search_dirs):
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

def run_combiner():
    parser = argparse.ArgumentParser(description="FASTQ Combiner (Improved)")
    parser.add_argument('csv_file', help='CSV mapping file')
    parser.add_argument('-o', '--output-dir', default='combined', help='Output directory')
    parser.add_argument('-d', '--search-dirs', nargs='+', help='Dirs to search for FASTQ')
    parser.add_argument('--validate-only', action='store_true', help='Validate R1/R2 read counts without combining')
    args = parser.parse_args()

    print(f"Running FASTQ Combiner Improved")
    print(f"CSV: {args.csv_file}")
    print(f"Output dir: {args.output_dir}")
    print(f"Search dirs: {args.search_dirs}")
    print(f"Validate only: {args.validate_only}")

    os.makedirs(args.output_dir, exist_ok=True)

    mapping = read_mapping_file(args.csv_file)
    print(f"Read mapping: {len(mapping)} targets")

    file_pairs = find_fastq_pairs(args.search_dirs)
    available_samples = list(file_pairs.keys())
    print(f"Discovered FASTQ pairs: {len(available_samples)} samples")

    results = []

    for target, sources in mapping.items():
        matched = []
        for s in sources:
            if s in file_pairs:
                matched.append(s)
            else:
                print(f"[WARNING] Source {s} not found for target {target}")
        if not matched:
            print(f"[WARNING] No valid sources found for target {target}")
            continue

        print(f"[{target}] Validating {len(matched)} source(s)")
        total_r1 = 0
        total_r2 = 0

        for src in matched:
            r1_file = file_pairs[src]['R1']
            r2_file = file_pairs[src]['R2']
            r1_count = count_reads_fastq(r1_file)
            r2_count = count_reads_fastq(r2_file)
            total_r1 += r1_count
            total_r2 += r2_count
            print(f"  {src}: R1={r1_count} R2={r2_count}")

        status = "PASS" if total_r1 == total_r2 else "FAIL"

        # If combining is enabled, produce combined FASTQs
        combined_r1_output = ""
        combined_r2_output = ""
        md5_r1 = ""
        md5_r2 = ""

        if not args.validate_only:
            clean_target = target.replace(" ", "_").replace("-", "_")
            combined_r1_output = os.path.join(args.output_dir, f"{clean_target}_S1_R1_001.fastq.gz")
            combined_r2_output = os.path.join(args.output_dir, f"{clean_target}_S1_R2_001.fastq.gz")
            print(f"[{target}] Combining R1 → {combined_r1_output}")
            print(f"[{target}] Combining R2 → {combined_r2_output}")
            r1_combined_reads = combine_fastq_files([file_pairs[s]['R1'] for s in matched], combined_r1_output)
            r2_combined_reads = combine_fastq_files([file_pairs[s]['R2'] for s in matched], combined_r2_output)
            md5_r1 = md5sum(combined_r1_output)
            md5_r2 = md5sum(combined_r2_output)
            print(f"[{target}] Combined R1 reads: {r1_combined_reads}, md5: {md5_r1}")
            print(f"[{target}] Combined R2 reads: {r2_combined_reads}, md5: {md5_r2}")

        results.append({
            "sample": target,
            "r1_count": total_r1,
            "r2_count": total_r2,
            "status": status,
            "combined_r1_output": combined_r1_output,
            "combined_r2_output": combined_r2_output,
            "md5_r1": md5_r1,
            "md5_r2": md5_r2
        })

    # Write HTML report
    generate_html_report(args.output_dir, results)

    print(f"Done. {len(results)} targets processed.")
