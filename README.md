# FASTQ Combiner
<img src="https://github.com/user-attachments/assets/c2ad7861-90f3-4448-b3a1-155245dd449e" alt="fastq_combiner_image" width="200" height="200">

A high-speed, memory-efficient tool for combining paired-end FASTQ files (supports both .fastq and .fastq.gz) for Cell Ranger and other NGS pipelines.

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](CONTRIBUTING.md)

---

## Features

- Streaming I/O for minimal RAM usage
- Supports both compressed (.fastq.gz) and uncompressed (.fastq) files
- Fuzzy sample name matching
- Parallel processing of samples
- Progress bar for combination process
- Logging with configurable verbosity
- Customizable buffer size
- Overwrite protection for output files (`--force` flag)
- Dry run mode for validation
- Outputs both HTML and CSV summary reports
- Custom R1/R2 file pattern support

![carbon](https://github.com/user-attachments/assets/115748c7-09b4-4dcc-81a5-64885cf9c9bc)

---

## Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/asomohammed/fastq-combiner.git
cd fastq-combiner

# Make executable
chmod +x fastq_combiner.py

# Ready to use!
python3 fastq_combiner.py --help
```

### Basic Usage

1. **Create a mapping CSV file:**
```csv
target_sample,source_file1,source_file2,source_file3
Patient001,Sample1_L001,Sample1_L002,Sample1_L003
Patient002,Sample2_L001,Sample2_L002,Sample2_L003
Control_Group,Control1,Control2,Negative_Control
```

2. **Run Combiner:**
```bash
python3 fastq_combiner.py mapping.csv
```

3. **Get Cell Ranger ready files:**
```
combined/
├── Patient001_S1_R1_001.fastq.gz
├── Patient001_S1_R2_001.fastq.gz
├── Patient002_S1_R1_001.fastq.gz
├── Patient002_S1_R2_001.fastq.gz
├── Control_Group_S1_R1_001.fastq.gz
├── Control_Group_S1_R2_001.fastq.gz
└── combination_report.html
```

---

## Detailed Usage

### Command Line Options

```bash
python3 fastq_combiner.py mapping.csv [options]

Options:
  -o, --output-dir DIR    Output directory (default: combined)
  -d, --search-dirs DIR   Directories to search for FASTQ files
  -h, --help             Show help message
```

### With Search Directories
```bash
python3 fastq_combiner.py mapping.csv -d /data/run1 /data/run2
```

### Custom Output Directory
```bash
python3 fastq_combiner.py mapping.csv -o cellranger_input
```

### Custom Buffer Size
```bash
python3 fastq_combiner.py mapping.csv --buffer-size 16777216  # 16MB
```

### Dry Run (Validate Only)
```bash
python3 fastq_combiner.py mapping.csv --dry-run
```

### Overwrite Output Files
```bash
python3 fastq_combiner.py mapping.csv --force
```

### Custom R1/R2 Patterns
```bash
python3 fastq_combiner.py mapping.csv --r1-pattern "*_R1_*.fq.gz" --r2-pattern "*_R2_*.fq.gz"
```

### Use as a CLI Tool (after pip install .)
```bash
pip install .
fastq-combiner mapping.csv -d /data/run1 /data/run2
```

### Use with a YAML Config File
Example `config.yaml`:
```yaml
csv_file: mapping.csv
output_dir: combined
search_dirs:
  - /data/run1
  - /data/run2
buffer_size: 16777216
threads: 8
force: true
r1_pattern:
  - "*_R1_*.fastq.gz"
r2_pattern:
  - "*_R2_*.fastq.gz"
```
Run with:
```bash
python3 fastq_combiner.py --config config.yaml
```

### Docker Usage
Build the image:
```bash
docker build -t fastq-combiner .
```
Run the tool (mount your data):
```bash
docker run --rm -v /path/to/data:/data fastq-combiner mapping.csv -d /data
```

### Diagnostics and Profiling
- Run environment diagnostics:
  ```bash
  python3 fastq_combiner.py --diagnostics
  ```
- Enable performance profiling:
  ```bash
  python3 fastq_combiner.py mapping.csv --profile
  ```

### CSV Format

The tool accepts flexible CSV formats:

#### Option 1: With header
```csv
target_sample,source_file1,source_file2,source_file3
Patient_001,Sample1_L001,Sample1_L002,Sample1_L003
Patient_002,Sample2_L001,Sample2_L002,Sample2_L003
```

#### Option 2: Without header
```csv
Patient_001,Sample1_L001,Sample1_L002,Sample1_L003
Patient_002,Sample2_L001,Sample2_L002,Sample2_L003
```

#### Option 3: Variable number of sources
```csv
target_sample,source_file1,source_file2,source_file3,source_file4
Patient_001,Sample1_L001,Sample1_L002,Sample1_L003,Sample1_L004
Patient_002,Sample2_L001,Sample2_L002
Control,Control_Sample
```

---

## Fuzzy Matching

The tool automatically handles common naming variations:

| CSV Entry   | Actual File                              | Match Type       |
|-------------|------------------------------------------|------------------|
| `sample1`   | `Sample1_S1_L001_R1_001.fastq.gz`        | Case insensitive |
| `Patient-1` | `Patient_1_S1_L001_R1_001.fastq.gz`      | Dash/underscore  |
| `ctrl`      | `Control_Sample_S1_L001_R1_001.fastq.gz` | Partial match    |
| `sample_a`  | `SampleA_S1_L001_R1_001.fastq.gz`        | Format variation |

## Supported File Formats

The tool automatically detects and handles various FASTQ naming conventions:

### Illumina Standard
```
Sample1_S1_L001_R1_001.fastq.gz
Sample1_S1_L001_R2_001.fastq.gz
```

### Multi-lane Illumina
```
Sample1_S1_L001_R1_001.fastq.gz
Sample1_S1_L002_R1_001.fastq.gz
Sample1_S1_L003_R1_001.fastq.gz
Sample1_S1_L004_R1_001.fastq.gz
```

### Simple naming
```
Sample1_R1.fastq.gz
Sample1_R2.fastq.gz
```

### Numbered pairs
```
Sample1_1.fastq.gz
Sample1_2.fastq.gz
```

### Dot notation
```
Sample1.R1.fastq.gz
Sample1.R2.fastq.gz
```

## Cell Ranger Integration

Output files follow Illumina's standard naming convention as required by Cell Ranger:

```bash
# Input files (any format)
Sample1_L001_R1_001.fastq.gz
Sample1_L002_R1_001.fastq.gz
Sample1_L003_R1_001.fastq.gz

# Output files (Cell Ranger compatible)
Sample1_S1_R1_001.fastq.gz
Sample1_S1_R2_001.fastq.gz

# Ready for Cell Ranger
cellranger count --id=sample1 --fastqs=combined/ --sample=Sample1
```

---

## Performance

### Speed Benchmarks
- **Sequential processing**: ~5-10 minutes per GB
- **Parallel processing**: ~2-3 seconds per GB

### Memory Usage
- **Scales with file size**: ~100MB RAM per GB of FASTQ data

### File Size Limits
- **Large file size support**: Tested with files up to 50GB
- **Handles thousands of files**: Recursive scanning

---

## Output Reports

The tool generates detailed HTML reports showing:

- **Combination statistics**: Read counts, file sizes
- **Fuzzy matches applied**: Shows what corrections were made
- **Processing summary**: Success/failure rates
- **File discovery details**: What files were found and used
- **Warnings and errors**: Missing files, mismatched read counts

---

## Advanced Usage

### Complex Directory Structures
```bash
# Handle complex nested directories
python3 fastq_combiner.py mapping.csv -d \
  /data/2023/run1/ \
  /data/2023/run2/ \
  /backup/sequencing_data/ \
  /external_drive/fastq_files/
```

### Batch Processing
```bash
# Process multiple mapping files
for mapping in *.csv; do
    echo "Processing $mapping..."
    python3 fastq_combiner.py "$mapping" -d /data/ -o "output_${mapping%.csv}"
done
```

### Integration with Workflows
```bash
#!/bin/bash
# Combine FASTQ files then run Cell Ranger
python3 fastq_combiner.py samples.csv -d /data/sequencing/ -o cellranger_input/

# Run Cell Ranger on combined files
for sample in $(cut -d',' -f1 samples.csv | tail -n +2); do
    cellranger count \
        --id="${sample}" \
        --fastqs=cellranger_input/ \
        --sample="${sample}" \
        --transcriptome=/path/to/reference/
done
```

---

## Troubleshooting & FAQ

**Q: I get `ModuleNotFoundError: No module named 'tqdm'` or `No module named 'yaml'`.**
- A: Run `pip install -r requirements.txt` to install all dependencies.

**Q: The script says 'No FASTQ files found!' but my files exist.**
- A: Make sure you use the `-d` flag to specify the correct search directory.

**Q: Permission denied or unreadable file errors?**
- A: Check file permissions. The script will skip files it cannot read.

**Q: How do I use a config file?**
- A: See the YAML example above and use `--config config.yaml`.

**Q: How do I run in Docker?**
- A: See the Docker usage example above. Mount your data with `-v`.

**Q: How do I check my environment?**
- A: Use `--diagnostics` to print Python version, dependencies, disk space, and more.

---

For more details, see the script's help:
```bash
python3 fastq_combiner.py --help
```

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Author

**Aso Omer Mohammed**  
3DBM, Neurosurgery, University Hospital Freiburg  
GitHub: [@asomohammed](https://github.com/asomohammed)  
Email: aso.mohammed@uniklinik-freiburg.de
