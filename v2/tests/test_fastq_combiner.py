import gzip
import os
from v2.fastq_combiner.utils import count_reads_fastq
import subprocess
import time
import stat

def generate_synthetic_fastq(path, num_reads=10):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, "wt") as f:
        for i in range(num_reads):
            f.write(f"@SEQ_ID_{i}\n")
            f.write("ACGTACGTACGT\n")
            f.write("+\n")
            f.write("FFFFFFFFFFFF\n")

def test_count_reads_fastq(tmp_path):
    # Generate R1 and R2 synthetic FASTQ
    r1_path = tmp_path / "TestSample_R1.fastq.gz"
    r2_path = tmp_path / "TestSample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=15)
    generate_synthetic_fastq(r2_path, num_reads=15)

    r1_count = count_reads_fastq(str(r1_path))
    r2_count = count_reads_fastq(str(r2_path))

    print(f"R1 count: {r1_count}, R2 count: {r2_count}")

    # Assert they match expected counts
    assert r1_count == 15
    assert r2_count == 15

def test_empty_fastq(tmp_path):
    # Generate empty R1 and R2 FASTQ files
    r1_path = tmp_path / "EmptySample_R1.fastq.gz"
    r2_path = tmp_path / "EmptySample_R2.fastq.gz"
    with gzip.open(r1_path, "wt") as f:
        pass
    with gzip.open(r2_path, "wt") as f:
        pass
    from v2.fastq_combiner.utils import count_reads_fastq
    r1_count = count_reads_fastq(str(r1_path))
    r2_count = count_reads_fastq(str(r2_path))
    assert r1_count == 0
    assert r2_count == 0

def test_corrupt_fastq(tmp_path):
    # Create a corrupt/truncated FASTQ file (incomplete record)
    r1_path = tmp_path / "CorruptSample_R1.fastq.gz"
    with gzip.open(r1_path, "wt") as f:
        f.write("@SEQ_ID_1\nACGTACGTACGT\n+\nFFFFFFFFFFFF\n")  # complete record
        f.write("@SEQ_ID_2\nACGTACGTACGT\n+\n")  # incomplete record (missing quality)
    from v2.fastq_combiner.utils import count_reads_fastq
    r1_count = count_reads_fastq(str(r1_path))
    # Should count only the complete record
    assert r1_count == 1

def test_missing_r2(tmp_path):
    # Create an R1 file without a corresponding R2 file
    r1_path = tmp_path / "LoneSample_R1.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=10)
    # Don't create the R2 file
    from v2.fastq_combiner.utils import count_reads_fastq
    r1_count = count_reads_fastq(str(r1_path))
    assert r1_count == 10
    # The file discovery logic should not include this R1 file since it has no R2 pair
    # This test verifies that the R1 file exists and is readable, but would be excluded
    # from the final file pairs due to missing R2

def test_mismatched_read_counts(tmp_path):
    # Create R1 and R2 files with different numbers of reads
    r1_path = tmp_path / "MismatchSample_R1.fastq.gz"
    r2_path = tmp_path / "MismatchSample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=15)
    generate_synthetic_fastq(r2_path, num_reads=10)  # Different count
    from v2.fastq_combiner.utils import count_reads_fastq
    r1_count = count_reads_fastq(str(r1_path))
    r2_count = count_reads_fastq(str(r2_path))
    assert r1_count == 15
    assert r2_count == 10
    assert r1_count != r2_count  # Verify they are indeed mismatched

def test_dry_run_no_output(tmp_path):
    # Create synthetic FASTQ files
    r1_path = tmp_path / "DryRunSample_R1.fastq.gz"
    r2_path = tmp_path / "DryRunSample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=5)
    generate_synthetic_fastq(r2_path, num_reads=5)
    # Create mapping CSV
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"DryRunSample,{r1_path},{r2_path}\n")
    # Run the script with --dry-run
    output_dir = tmp_path / "dryrun_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--dry-run"
    ], capture_output=True, text=True)
    # No output files should be created
    assert not any(output_dir.glob("*.fastq.gz")), "No output files should be created in dry run mode"
    # Script should complete without error
    assert result.returncode == 0

def test_overwrite_protection(tmp_path):
    # Create synthetic FASTQ files
    r1_path = tmp_path / "OverwriteSample_R1.fastq.gz"
    r2_path = tmp_path / "OverwriteSample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=5)
    generate_synthetic_fastq(r2_path, num_reads=5)
    # Create mapping CSV
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"OverwriteSample,{r1_path},{r2_path}\n")
    # Run the script to generate output files
    output_dir = tmp_path / "overwrite_output"
    subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--force"
    ], check=True)
    r1_out = output_dir / "OverwriteSample_S1_R1_001.fastq.gz"
    r2_out = output_dir / "OverwriteSample_S1_R2_001.fastq.gz"
    assert r1_out.exists() and r2_out.exists()
    # Record modification times
    r1_mtime_before = r1_out.stat().st_mtime
    r2_mtime_before = r2_out.stat().st_mtime
    time.sleep(1)  # Ensure mtime would change if overwritten
    # Run the script again without --force
    subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path)
    ], check=True)
    # Modification times should not change
    assert r1_out.stat().st_mtime == r1_mtime_before
    assert r2_out.stat().st_mtime == r2_mtime_before

def test_permission_error(tmp_path):
    # Create a synthetic FASTQ file
    r1_path = tmp_path / "NoPermSample_R1.fastq.gz"
    r2_path = tmp_path / "NoPermSample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=5)
    generate_synthetic_fastq(r2_path, num_reads=5)
    # Remove read permissions from R1
    r1_path.chmod(0)
    # Create mapping CSV
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"NoPermSample,{r1_path},{r2_path}\n")
    # Run the script (should log an error and skip the file)
    output_dir = tmp_path / "noperm_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path)
    ], capture_output=True, text=True)
    # Restore permissions for cleanup
    r1_path.chmod(stat.S_IRUSR | stat.S_IWUSR)
    # Should log an error about unreadable file
    assert "Error" in result.stdout or "error" in result.stderr or "unreadable" in result.stdout.lower()
    # No output files should be created
    assert not any(output_dir.glob("*.fastq.gz")), "No output files should be created if input is unreadable"

def test_paired_end_deduplication(tmp_path):
    # Create synthetic FASTQ files with duplicate reads
    r1_path1 = tmp_path / "Sample1_R1.fastq.gz"
    r2_path1 = tmp_path / "Sample1_R2.fastq.gz"
    r1_path2 = tmp_path / "Sample2_R1.fastq.gz"
    r2_path2 = tmp_path / "Sample2_R2.fastq.gz"
    # Both samples have the same reads (simulate duplicates)
    generate_synthetic_fastq(r1_path1, num_reads=5)
    generate_synthetic_fastq(r2_path1, num_reads=5)
    generate_synthetic_fastq(r1_path2, num_reads=5)
    generate_synthetic_fastq(r2_path2, num_reads=5)
    # Create mapping CSV
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"DedupSample,{r1_path1},{r1_path2}\n")
    # Run the script with paired-end deduplication
    output_dir = tmp_path / "dedup_output"
    subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--deduplicate", "--paired-end-dedup", "--force"
    ], check=True)
    r1_out = output_dir / "DedupSample_S1_R1_001.fastq.gz"
    r2_out = output_dir / "DedupSample_S1_R2_001.fastq.gz"
    # Only 1 unique read pair should remain (since all reads are identical)
    from v2.fastq_combiner.utils import count_reads_fastq
    assert count_reads_fastq(str(r1_out)) == 1
    assert count_reads_fastq(str(r2_out)) == 1

def test_checksum_and_error_reporting(tmp_path):
    # Create synthetic FASTQ files
    r1_path = tmp_path / "CSVSamp_R1.fastq.gz"
    r2_path = tmp_path / "CSVSamp_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=3)
    generate_synthetic_fastq(r2_path, num_reads=3)
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"CSVSamp,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "csv_output"
    subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--force"
    ], check=True)
    # Check CSV summary for checksums and no errors
    csv_summary = output_dir / "combination_summary.csv"
    with open(csv_summary) as f:
        lines = f.readlines()
    assert "R1 Checksum" in lines[0]
    assert "R2 Checksum" in lines[0]
    assert "CSVSamp" in lines[1]
    # Allow ',0,' if it means no errors/skipped; only fail if skipped is True or error is non-empty
    columns = lines[1].split(',')
    skipped_col = columns[7].strip().lower()
    error_col = columns[8].strip().lower()
    assert skipped_col == 'false'
    assert error_col == ''

def test_quality_score_format_detection(tmp_path):
    # Create a FASTQ file with Sanger quality scores
    r1_path = tmp_path / "QualSample_R1.fastq"
    with open(r1_path, "w") as f:
        for i in range(5):
            f.write(f"@SEQ_ID_{i}\nACGTACGTACGT\n+\n!!!!!#####\n")  # ASCII 33-35
    r2_path = tmp_path / "QualSample_R2.fastq"
    with open(r2_path, "w") as f:
        for i in range(5):
            f.write(f"@SEQ_ID_{i}\nACGTACGTACGT\n+\n!!!!!#####\n")
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"QualSample,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "qual_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--validate", "--force"
    ], capture_output=True, text=True)
    # Should log detected quality format
    assert "Detected quality format" in result.stdout or "Detected quality format" in result.stderr

def test_paired_end_integrity_error(tmp_path):
    # Create R1 and R2 with mismatched reads
    r1_path = tmp_path / "PEInt_R1.fastq.gz"
    r2_path = tmp_path / "PEInt_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=5)
    generate_synthetic_fastq(r2_path, num_reads=3)
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"PEInt,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "peint_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path)
    ], capture_output=True, text=True)
    # Should log a paired-end mismatch warning or error
    assert "paired-end" in result.stdout.lower() or "paired-end" in result.stderr.lower()

def test_adapter_and_gc_analysis(tmp_path):
    # Create FASTQ with known adapter and GC content
    r1_path = tmp_path / "GCAdapter_R1.fastq"
    r2_path = tmp_path / "GCAdapter_R2.fastq"
    with open(r1_path, "w") as f:
        for i in range(5):
            f.write(f"@SEQ_ID_{i}\nACGTACGTACGTAGATCGGAAGAGC\n+\nFFFFFFFFFFFFFFFFFFFFF\n")  # Adapter
    with open(r2_path, "w") as f:
        for i in range(5):
            f.write(f"@SEQ_ID_{i}\nGCGCGCGCGCGC\n+\nFFFFFFFFFFFF\n")
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"GCAdapter,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "gc_adapter_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--adapter-check", "--gc-analysis", "--force"
    ], capture_output=True, text=True)
    # Should log adapter detection and GC content
    assert "adapter" in result.stdout.lower() or "adapter" in result.stderr.lower()
    assert "gc content" in result.stdout.lower() or "gc content" in result.stderr.lower()

def test_fuzzy_matching(tmp_path):
    # Create FASTQ files with slightly different sample names
    r1_path = tmp_path / "FuzzySample_R1.fastq.gz"
    r2_path = tmp_path / "FuzzySample_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=2)
    generate_synthetic_fastq(r2_path, num_reads=2)
    mapping_csv = tmp_path / "mapping.csv"
    # Intentionally typo in mapping
    with open(mapping_csv, "w") as f:
        f.write(f"FuzzySampel,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "fuzzy_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--force"
    ], capture_output=True, text=True)
    # Should log a fuzzy match
    assert "fuzzy match" in result.stdout.lower() or "fuzzy match" in result.stderr.lower()

def test_dry_run_all_features(tmp_path):
    # Create synthetic FASTQ files
    r1_path = tmp_path / "DryAll_R1.fastq.gz"
    r2_path = tmp_path / "DryAll_R2.fastq.gz"
    generate_synthetic_fastq(r1_path, num_reads=2)
    generate_synthetic_fastq(r2_path, num_reads=2)
    mapping_csv = tmp_path / "mapping.csv"
    with open(mapping_csv, "w") as f:
        f.write(f"DryAll,{r1_path},{r2_path}\n")
    output_dir = tmp_path / "dryall_output"
    result = subprocess.run([
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "--search-dirs", str(tmp_path), "--dry-run", "--deduplicate", "--paired-end-dedup", "--validate", "--adapter-check", "--gc-analysis"
    ], capture_output=True, text=True)
    # Should not create any output files
    assert not any(output_dir.glob("*.fastq.gz")), "No output files should be created in dry run mode"
    # Should log dry run
    assert "dry run" in result.stdout.lower() or "dry run" in result.stderr.lower()
