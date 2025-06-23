import gzip
import os
import random

bases = 'ACGT'
def write_fastq(path, n):
    if path.endswith('.gz'):
        f = gzip.open(path, 'wt')
    else:
        f = open(path, 'w')
    for i in range(n):
        seq = ''.join(random.choices(bases, k=100))
        qual = 'I' * 100
        f.write(f'@SEQ{i}\n{seq}\n+\n{qual}\n')
    f.close()

dirs = {
    'run1': 'Sample_A',
    'run2': 'Sample_A',
    'run3': 'Sample_A',
    'run4': 'Sample_A',
    'batch1/data': 'Sample_B',
    'batch2/data': 'Sample_B',
    'batch3/data': 'Sample_B',
    'local': 'Sample_C',
}

for d, s in dirs.items():
    os.makedirs(d, exist_ok=True)
    write_fastq(f'{d}/{s}_R1.fastq.gz', 10)
    write_fastq(f'{d}/{s}_R2.fastq.gz', 10) 