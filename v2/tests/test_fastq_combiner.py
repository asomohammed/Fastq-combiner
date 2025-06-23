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
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "-d", str(tmp_path), "--force"
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
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "-d", str(tmp_path)
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
        "python3", "fastq_combiner.py", str(mapping_csv), "-o", str(output_dir), "-d", str(tmp_path)
    ], capture_output=True, text=True)
    # Restore permissions for cleanup
    r1_path.chmod(stat.S_IRUSR | stat.S_IWUSR)
    # Should log an error about unreadable file
    assert "Error" in result.stdout or "error" in result.stderr or "unreadable" in result.stdout.lower()
    # No output files should be created
    assert not any(output_dir.glob("*.fastq.gz")), "No output files should be created if input is unreadable"
