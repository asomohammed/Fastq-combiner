import gzip

def count_reads_fastq(fastq_file: str) -> int:
    count = 0
    opener = gzip.open if fastq_file.endswith('.gz') else open
    with opener(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0 and line.startswith('@'):
                count += 1
    return count
