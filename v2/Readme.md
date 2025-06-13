## FASTQ Combiner v2
[![Python](https://img.shields.io/badge/Python-3.8%2B-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Build Status](https://img.shields.io/badge/status-Stable-brightgreen.svg)]()

---

## Overview

**FASTQ Combiner Pipeline is an easy-to-use tool for combining FASTQ sequencing files and validating the results.**

When sequencing large projects, samples are often split across multiple lanes or batches. This tool allows you to:
- Automatically find matching pairs of FASTQ files (R1 and R2)
- Combine FASTQ files from multiple sources into a single output FASTQ pair for each sample
- Validate that your R1 and R2 files match (important for quality control)
- Generate a clean HTML report summarizing the results
- Use a “validate-only” mode to quickly check your files before combining them

---

## Folder Structure
```text
fastq_combiner/v2
├── __main__.py         # CLI entry point
├── core.py             # Main pipeline logic
├── report.py           # Bootstrap 5 HTML report generator
├── utils.py            # Utility functions (R1/R2 read counting, md5sum, combine_fastq_files)
tests/
└── test_fastq_combiner.py   # Unit test scaffold (synthetic FASTQ test)
requirements.txt        # Dependencies
README.md               # This file
.gitignore              # Git ignore rules
.github/workflows/tests.yml # GitHub Actions automated testing workflow
```

---

## Installation

1. Download the pipeline

Clone the repository or download it as a ZIP and extract it:
```bash
git clone https://github.com/YOUR_USERNAME/fastq_combiner.git
cd fastq_combiner
```
2. Install required Python packages

You need Python 3.8 or higher.

Install the required Python packages:
```bash
pip install -r requirements.txt
```

---

## Working Concept
- You provide a CSV file that maps your target samples to one or more source sample names.
- The tool scans your provided directories and automatically finds matching FASTQ files.
- It can either:
  - Just validate that R1 and R2 files match (validate-only mode)
  - Or fully combine the files into a new FASTQ pair (combine mode)
  - It produces a clean HTML report to summarize what it did.

---

## Usage

### Validate-only mode

This will not create any new FASTQ files. It only checks that your input files are valid.

```bash
python3 -m fastq_combiner mapping.csv -o combined_output -d /path/to/fastq_dir1 /path/to/fastq_dir2 --validate-only
```

What it does:
- Scans the FASTQ directories you provide
- Looks for matching pairs of R1 and R2 FASTQ files
- For each target sample in your CSV, finds the matching sources
- Counts the reads in R1 and R2 and checks that they match
- Generates an HTML report at:
```bash
combined_output/combination_report.html
```

### Full combine mode

```bash
python3 -m fastq_combiner mapping.csv -o combined_output -d /path/to/fastq_dir1 /path/to/fastq_dir2
```

What it does:
- Same steps as validate-only
- PLUS combines the matching FASTQ files into a new R1 and R2 FASTQ pair per target sample
- Outputs the new combined FASTQ files to your combined_output/ directory
- Computes MD5 checksums for reproducibility

---


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Citation

If you use this tool in your work, please consider citing:

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
