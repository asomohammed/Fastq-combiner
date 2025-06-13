## Overview

FASTQ Combiner is a robust, parallelized, and resumable FASTQ file merging pipeline. It enables flexible combination of multiple source FASTQ pairs (R1/R2) into target sample FASTQ pairs, guided by a user-provided CSV mapping.

Features:
- Automatic FASTQ pair discovery (supports diverse FASTQ naming conventions)
- Fuzzy matching of sample names
- Parallel processing with thread pools
- Resume support (skip already existing outputs)
- Detailed integrity checks (MD5 checksums, read count consistency)
- Generates HTML, CSV, and JSON reports

Designed for managing large-scale sequencing projects where samples may be split across multiple lanes or batches.

---

## Features
- Automatic FASTQ pair discovery
- Fuzzy matching of sample names (optional fallback to difflib)
- Parallel processing with configurable number of jobs
- Resume support
- Robust logging and reporting
- Generates HTML, CSV, and JSON summary outputs
- Supports gzipped and uncompressed FASTQ
- Efficient I/O with buffer tuning

---

## Installation
```bash
git clone https://github.com/YOUR_USERNAME/fastq_combiner.git
cd fastq_combiner
pip install -r requirements.txt
pip install tqdm psutil rapidfuzz
```

---

## Usage
```bash
python3 fastq_combiner.py mapping.csv -o combined_output_dir -d /path/to/fastq_dir1 /path/to/fastq_dir2 --resume --log INFO
```

### Required argument
- mapping.csv: CSV file defining the mapping of target sample name to source sample names.

### Example:

TargetSample1, SourceSample1A, SourceSample1B, SourceSample1C
TargetSample2, SourceSample2A, SourceSample2B

### Common options

|Option	          |    Description	                                      |     Default         |
|-----------------|-------------------------------------------------------|---------------------|
|-o, --output-dir | Output directory for combined FASTQ files and reports	| combined            |
|-d, --search-dirs|	 Directories to recursively search for FASTQ files	  |  Current directory  |
|--resume	        |  Skip targets where output files already exist	      |  Off                |
|--dry-run	      |  Simulate run without producing combined files	      |  Off                |
|--buffer-size    |  Buffer size in bytes for combining files	            |  8 MB               |
|--compresslevel  |  Gzip compression level for outputs	                  |  6                  |
|-j, --jobs	      |  Number of parallel threads to use	                  |  All available CPUs |
|--log-file	      |  Log file path (optional)	                            |  None               |
|--log	          |  Log level (INFO, DEBUG, WARNING, etc.)	              |  INFO               |

---

## Outputs

For each target sample, the script produces:
- Combined *_R1_001.fastq.gz and *_R2_001.fastq.gz FASTQ files
- summary.json — full processing metadata
- summary.csv — tabular summary
- combination_report.html — interactive HTML report with read counts, status, and integrity checks

⸻

## Reporting

The HTML report provides:
- Total number of targets processed
- Total number of reads combined
- Per-target status (Success/Failed)
- Read counts per target
- Fuzzy matches used (if applicable)
- Links to summary CSV/JSON files

---

## License

MIT License

---

## Citation

If you use this tool in your work, please cite it appropriately:

FASTQ Combiner v2.0, https://github.com/asomohammed/fastq_combiner

---

## Authors

**Aso Omer Mohammed**  
3DBM, Neurosurgery, University Hospital Freiburg  
GitHub: [@asomohammed](https://github.com/asomohammed)  
Email: aso.mohammed@uniklinik-freiburg.de

**Dr. Ing. Kevin Joseph**  
Neurosurgery, University Hospital Freiburg  
GitHub: [@kevinj24fr](https://github.com/kevinj24fr)  
Email: kevin.joseph@uniklinik-freiburg.de
