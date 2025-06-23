# ⚡ FASTQ File Combiner - STREAMING OPTIMIZED

**High-speed streaming I/O with minimal RAM usage for combining paired-end FASTQ files.**
Optimized for Cell Ranger compatibility and large-scale sequencing data.

## 🚀 Key Features

### **Performance & Scalability**
- **Streaming I/O**: Minimal RAM usage, handles files of any size
- **Auto-optimization**: Automatically tunes buffer sizes based on system resources
- **Memory profiling**: Track RAM usage during processing
- **Checkpointing**: Resume interrupted runs from where they left off
- **SSD/HDD detection**: Optimizes I/O patterns for different storage types
- **Multi-threading**: Parallel processing for maximum speed

### **Data Quality & Validation**
- **Quality score validation**: Detects corrupted quality scores
- **Read length consistency**: Verifies R1/R2 have matching lengths
- **Adapter detection**: Identifies common adapter contamination
- **GC content analysis**: Basic sequence quality metrics
- **Sample barcode extraction**: Parse and validate sample barcodes

### **Monitoring & Safety**
- **Real-time monitoring**: Live updates of processing speed/throughput
- **Disk space monitoring**: Warns before running out of space
- **Backup creation**: Auto-backup original files before processing
- **Retry logic**: Automatically retries failed operations
- **Corruption detection**: Validates file integrity before processing

### **User Experience**
- **Progress estimation**: Shows ETA for large datasets
- **Interactive HTML reports**: Search, filter, and visualize results
- **Configuration files**: YAML-based configuration for complex workflows
- **System diagnostics**: Comprehensive system information
- **Performance profiling**: Detailed performance analysis

## 📦 Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/Fastq-combiner.git
cd Fastq-combiner

# Install dependencies
pip install -r requirements.txt

# Or install directly
pip install tqdm pyyaml psutil
```

## 🎯 Quick Start

### Basic Usage
```bash
python3 fastq_combiner.py mapping.csv -o combined_output
```

### Enhanced Usage with All Features
```bash
python3 fastq_combiner.py mapping.csv \
  --config example_config.yaml \
  --auto-optimize \
  --memory-profile \
  --validate \
  --real-time-monitor \
  --checkpoint
```

## 📋 Input Format

### CSV Mapping File
```csv
target_sample,source_file1,source_file2,source_file3
Sample_A,run1/Sample_A,run2/Sample_A,run3/Sample_A
Sample_B,batch1/Sample_B,batch2/Sample_B,batch3/Sample_B
Sample_C,seq1/Sample_C,seq2/Sample_C,local/Sample_C
```

### YAML Configuration
```yaml
# example_config.yaml
output: "combined_output"
threads: 8
auto_optimize: true
validate: true
monitor_disk: true
create_backups: true
```

## 🔧 Advanced Features

### Performance Optimization
```bash
# Auto-optimize based on system resources
--auto-optimize

# Memory profiling
--memory-profile

# Performance profiling
--profile

# Custom buffer size
--buffer-size 33554432  # 32MB
```

### Data Validation
```bash
# Validate FASTQ quality
--validate

# Check sample barcodes
--check-barcodes

# GC content analysis
--gc-analysis

# Adapter detection
--adapter-check
```

### Monitoring & Safety
```bash
# Real-time monitoring
--real-time-monitor

# Disk space monitoring
--monitor-disk

# Create backups
--create-backups

# Retry failed operations
--retry-failed
```

### Checkpointing
```bash
# Enable checkpointing
--checkpoint

# Resume interrupted run
python3 fastq_combiner.py mapping.csv --checkpoint
```

## 📊 Output

### Files Generated
- `{sample}_S1_R1_001.fastq.gz` - Combined R1 reads
- `{sample}_S1_R2_001.fastq.gz` - Combined R2 reads
- `combination_summary.csv` - Processing summary
- `combination_report.html` - Interactive HTML report

### HTML Report Features
- **Interactive search/filter** for large datasets
- **Success/failure charts** with Chart.js
- **Per-sample details** with collapsible sections
- **Validation results** and quality metrics
- **System metadata** and performance stats

## 🐳 Docker Support

```bash
# Build image
docker build -t fastq-combiner .

# Run with volume mount
docker run -v $(pwd):/data fastq-combiner mapping.csv -o /data/output
```

## 🔍 System Diagnostics

```bash
# Print system information
python3 fastq_combiner.py --diagnostics
```

Output includes:
- Python version and platform
- CPU cores and memory
- Disk space and storage type
- Dependency status

## 📈 Performance Tips

1. **Use SSD storage** for best performance
2. **Enable auto-optimization** for automatic tuning
3. **Monitor memory usage** with `--memory-profile`
4. **Use checkpointing** for large datasets
5. **Enable real-time monitoring** for progress tracking

## 🛠️ Troubleshooting

### Common Issues

**"No FASTQ files found"**
- Check search directories with `--search-dirs`
- Verify file patterns with `--r1-patterns` and `--r2-patterns`

**"Low disk space"**
- Use `--monitor-disk` to check space requirements
- Clean up temporary files

**"Memory issues"**
- Use `--memory-profile` to monitor usage
- Reduce `--buffer-size` or `--threads`

**"Validation warnings"**
- Use `--validate` to check file quality
- Review warnings in HTML report

### Performance Optimization

```bash
# For large datasets
--auto-optimize --memory-profile --checkpoint --real-time-monitor

# For validation-heavy workflows
--validate --check-barcodes --gc-analysis --adapter-check

# For safety-critical operations
--create-backups --monitor-disk --retry-failed
```

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new features
4. Submit a pull request

See [CONTRIBUTING.md](CONTRIBUTING.md) for details.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Built for the bioinformatics community
- Optimized for Cell Ranger compatibility
- Inspired by the need for efficient large-scale data processing
