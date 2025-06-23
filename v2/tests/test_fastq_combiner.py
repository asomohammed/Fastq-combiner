import gzip
import os
from v2.fastq_combiner.utils import count_reads_fastq

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
