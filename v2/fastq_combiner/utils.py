import gzip
import hashlib
import shutil
import os

def count_reads_fastq(fastq_file: str) -> int:
    count = 0
    opener = gzip.open if fastq_file.endswith('.gz') else open
    with opener(fastq_file, 'rt') as f:
        while True:
            # Read all four lines of a FASTQ record
            header = f.readline()
            if not header:
                break
            seq = f.readline()
            plus = f.readline()
            qual = f.readline()
            # Only count if all four lines are present and valid
            if (header.startswith('@') and seq and plus and qual):
                count += 1
    return count

def md5sum(filename: str, blocksize: int = 2**20) -> str:
    m = hashlib.md5()
    with open(filename, "rb") as f:
        for block in iter(lambda: f.read(blocksize), b""):
            m.update(block)
    return m.hexdigest()

def get_decompressed_size(filename: str) -> int:
    """Get actual decompressed size of a gzipped file"""
    if not filename.endswith('.gz'):
        return os.path.getsize(filename)
    
    with gzip.open(filename, 'rb') as f:
        # Read in chunks to avoid memory issues
        size = 0
        while True:
            chunk = f.read(1024 * 1024)  # 1MB chunks
            if not chunk:
                break
            size += len(chunk)
    return size

def combine_fastq_files(source_files: list, output_file: str, buffer_size: int = 8 * 1024 * 1024, compresslevel: int = 6) -> int:
    total_reads = 0
    
    opener = gzip.open if output_file.endswith('.gz') else open
    with opener(output_file, 'wb', compresslevel=compresslevel) if output_file.endswith('.gz') \
         else open(output_file, 'wb') as out_f:
        
        for i, src in enumerate(source_files, 1):
            file_size = os.path.getsize(src)
            file_type = "compressed" if src.endswith('.gz') else "uncompressed"
            print(f"    [{i}/{len(source_files)}] Processing {os.path.basename(src)} ({file_size / (1024**3):.2f} GB {file_type})")
            
            src_opener = gzip.open if src.endswith('.gz') else open
            with src_opener(src, 'rb') as in_f:
                bytes_processed = 0
                while True:
                    chunk = in_f.read(buffer_size)
                    if not chunk:
                        break
                    out_f.write(chunk)
                    bytes_processed += len(chunk)
                    
                    # Show progress every 200MB - ONLY data processed, NO percentage
                    if bytes_processed % (200 * 1024 * 1024) < len(chunk):
                        print(f"      Processed: {bytes_processed / (1024**3):.2f} GB")
            
            print(f"      ✓ Complete: {os.path.basename(src)} - Total processed: {bytes_processed / (1024**3):.2f} GB")
            total_reads += count_reads_fastq(src)
    
    return total_reads
