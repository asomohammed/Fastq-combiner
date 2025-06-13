import gzip
import hashlib
import shutil

def count_reads_fastq(fastq_file: str) -> int:
    count = 0
    opener = gzip.open if fastq_file.endswith('.gz') else open
    with opener(fastq_file, 'rt') as f:
        for i, line in enumerate(f):
            if i % 4 == 0 and line.startswith('@'):
                count += 1
    return count

def md5sum(filename: str, blocksize: int = 2**20) -> str:
    m = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            m.update(block)
    return m.hexdigest()

def combine_fastq_files(source_files: list, output_file: str, buffer_size: int = 8 * 1024 * 1024, compresslevel: int = 6) -> int:
    total_reads = 0
    opener = gzip.open if output_file.endswith('.gz') else open
    with opener(output_file, 'wb', compresslevel=compresslevel) if output_file.endswith('.gz') \
         else open(output_file, 'wb') as out_f:
        for src in source_files:
            src_opener = gzip.open if src.endswith('.gz') else open
            with src_opener(src, 'rb') as in_f:
                shutil.copyfileobj(in_f, out_f, length=buffer_size)
            total_reads += count_reads_fastq(src)
    return total_reads
