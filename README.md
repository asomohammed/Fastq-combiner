# FASTQ Combiner

<img src="https://github.com/user-attachments/assets/c2ad7861-90f3-4448-b3a1-155245dd449e" alt="fastq_combiner_image" width="200" height="200">


A fast, intelligent FASTQ file combiner that generates Cell Ranger compatible outputs with fuzzy matching for sample names. Perfect for combining multi-lane sequencing data or merging samples across different sequencing runs.

[![Python 3.6+](https://img.shields.io/badge/python-3.6+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Contributions Welcome](https://img.shields.io/badge/contributions-welcome-brightgreen.svg)](CONTRIBUTING.md)

---

## Features

- Parallel processing with multi-threading
- Automatic Illumina-standard naming (`{sample}_S1_R1_001.fastq.gz`)
- Handles typos and naming variations automatically
- Recursive file scanning with intelligent pattern matching
-  HTML reports with combination statistics
- Seamlessly combines files from different lanes
- Supports various FASTQ naming conventions

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
python3 fastq_combiner.py mapping.csv -d /path/to/fastq/files
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

---

## Usage Example

#### Combine multi-lane data
```bash
# Combine lanes from a single run
python3 fastq_combiner.py samples.csv -d /data/NovaSeq_Run_001/
```

#### Combine across multiple runs
```bash
# Search multiple directories
python3 fastq_combiner.py samples.csv -d /data/run1/ /data/run2/ /data/run3/
```

#### Custom output directory
```bash
# Output to cellranger_input folder
python3 fastq_combiner.py samples.csv -o cellranger_input -d /data/
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

## Troubleshooting

### Common Issues

#### Files not found
```bash
# Check file discovery
python3 fastq_combiner.py mapping.csv -d /data/ 2>&1 | grep "Found"
```

#### Permission errors
```bash
# Fix permissions
chmod -R 755 /path/to/fastq/files/
```

#### Memory issues with large files
```bash
# Process smaller batches
split -l 10 large_mapping.csv batch_
for batch in batch_*; do
    python3 fastq_combiner.py "$batch" -d /data/
done
```

### Debug Mode
```bash
# Verbose output for debugging
python3 fastq_combiner.py mapping.csv -d /data/ 2>&1 | tee debug.log
```

### Development Setup
```bash
git clone https://github.com/yourusername/fastq-combiner.git
cd fastq-combiner

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/

# Format code
black fastq_combiner.py
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
