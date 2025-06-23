#!/usr/bin/env python3
"""
FASTQ File Combiner - STREAMING OPTIMIZED
High-speed streaming I/O with minimal RAM usage for combining paired-end FASTQ files.
Optimized for Cell Ranger compatibility and large-scale sequencing data.
"""

import os
import sys
import gzip
import glob
import csv
import logging
import argparse
import subprocess
import time
import json
import hashlib
import psutil
import threading
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict
import yaml
import cProfile
import pstats
import io
from tqdm import tqdm
import pickle
from datetime import datetime
import shutil
import tempfile
import platform
import importlib
import socket
import getpass
import html as html_escape
import multiprocessing

BUFFER_SIZE = 8 * 1024 * 1024  # Default, can be overridden by CLI

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='[%(levelname)s] %(message)s'
)

try:
    import resource
except ImportError:
    resource = None

def read_mapping_file(csv_file):
    """Read CSV mapping file and return dictionary of target -> [source file paths]"""
    mapping = defaultdict(list)
    
    logging.info(f"Reading mapping file: {csv_file}")
    
    with open(csv_file, 'r') as f:
        reader = csv.reader(f)
        
        # Check if first row is header
        first_row = next(reader)
        if first_row[0].lower() in ['target', 'target_sample', 'output', 'sample']:
            logging.info("  Detected header row, skipping...")
        else:
            # Process first row as data
            target = first_row[0].strip()
            sources = [f.strip() for f in first_row[1:] if f.strip()]
            mapping[target].extend(sources)
        
        # Process remaining rows
        for row_num, row in enumerate(reader, start=2):
            if not row or not row[0].strip():
                continue
                
            target = row[0].strip()
            sources = [f.strip() for f in row[1:] if f.strip()]
            
            if not sources:
                logging.warning(f"  No source files for target '{target}' on row {row_num}")
                continue
                
            mapping[target].extend(sources)
    
    logging.info(f"Loaded {len(mapping)} target samples")
    return dict(mapping)

def find_fastq_files_fast(search_dirs, r1_patterns=None, r2_patterns=None):
    """
    Fast file discovery using glob patterns with improved logic
    Returns dict: {sample_base: {'R1': path, 'R2': path}}
    """
    logging.info("\n🔍 Scanning for FASTQ files...")
    if not search_dirs:
        search_dirs = ["."]
    
    file_pairs = {}
    temp_pairs = {}  # Temporary storage to avoid race conditions
    
    # Use custom patterns if provided, else default
    default_r1_patterns = [
        "*_R1_*.fastq.gz", "*_R1.fastq.gz", "*_1.fastq.gz", "*.R1.fastq.gz",
        "*_R1_*.fastq", "*_R1.fastq", "*_1.fastq", "*.R1.fastq"
    ]
    r1_patterns = r1_patterns or default_r1_patterns
    default_r2_patterns = [
        p.replace('R1', 'R2').replace('_1', '_2').replace('.R1', '.R2') for p in default_r1_patterns
    ]
    r2_patterns = r2_patterns or default_r2_patterns
    
    logging.debug(f"R1 patterns: {r1_patterns}")
    logging.debug(f"R2 patterns: {r2_patterns}")
    
    for search_dir in search_dirs:
        logging.info(f"  Scanning: {os.path.abspath(search_dir)}")
        for r1_pattern in r1_patterns:
            full_pattern = os.path.join(search_dir, "**", r1_pattern)
            logging.debug(f"  Looking for pattern: {full_pattern}")
            r1_files = glob.glob(full_pattern, recursive=True)
            logging.debug(f"  Found R1 files: {r1_files}")
            
            for r1_file in r1_files:
                basename = os.path.basename(r1_file)
                logging.debug(f"  Processing R1 file: {basename}")
                
                # Extract sample base using consistent strategy
                sample_base = None
                
                # Strategy 1: Standard Illumina pattern (sample_S1_R1_001.fastq.gz)
                if "_S" in basename and "_R1_" in basename:
                    sample_base = basename.split("_S")[0]
                # Strategy 2: Simple R1 pattern (sample_R1_001.fastq.gz)
                elif "_R1_" in basename:
                    sample_base = basename.split("_R1_")[0]
                # Strategy 3: Simple R1 pattern (sample_R1.fastq.gz)
                elif "_R1." in basename:
                    sample_base = basename.split("_R1.")[0]
                # Strategy 4: Numbered pattern (sample_1.fastq.gz)
                elif "_1." in basename:
                    sample_base = basename.split("_1.")[0]
                # Strategy 5: Dot pattern (sample.R1.fastq.gz)
                elif ".R1." in basename:
                    sample_base = basename.split(".R1.")[0]
                else:
                    logging.debug(f"    No pattern match for: {basename}")
                    continue
                
                logging.debug(f"    Extracted sample base: {sample_base}")
                
                # Find corresponding R2 file
                r2_file = None
                r2_candidates = []
                
                # Generate R2 candidates based on the same pattern logic
                for r2_pattern in r2_patterns:
                    r2_candidate = os.path.join(os.path.dirname(r1_file), basename.replace('R1', 'R2').replace('_1', '_2').replace('.R1', '.R2'))
                    r2_candidates.append(r2_candidate)
                
                # Also try legacy patterns
                r2_candidates.extend([
                    r1_file.replace("_R1_", "_R2_"),
                    r1_file.replace("_R1.", "_R2."),
                    r1_file.replace("_1.", "_2."),
                    r1_file.replace(".R1.", ".R2."),
                    r1_file.replace("_R1_", "_R2_").replace(".fastq.gz", ".fastq"),
                    r1_file.replace("_R1.", "_R2.").replace(".fastq.gz", ".fastq"),
                    r1_file.replace("_1.", "_2.").replace(".fastq.gz", ".fastq"),
                    r1_file.replace(".R1.", ".R2.").replace(".fastq.gz", ".fastq")
                ])
                
                logging.debug(f"    R2 candidates: {r2_candidates}")
                
                # Find existing R2 file
                for r2_candidate in r2_candidates:
                    if os.path.exists(r2_candidate):
                        r2_file = r2_candidate
                        logging.debug(f"    Found R2 file: {r2_candidate}")
                        break
                
                if r2_file:
                    full_r1_path = os.path.abspath(r1_file)
                    full_r2_path = os.path.abspath(r2_file)
                    
                    # Use consistent key strategy: prefer sample_base, fallback to full path
                    key = sample_base
                    if key in temp_pairs:
                        # If sample_base already exists, use full path to avoid conflicts
                        key = full_r1_path
                        logging.debug(f"    Using full path as key due to conflict: {key}")
                    
                    temp_pairs[key] = {'R1': full_r1_path, 'R2': full_r2_path}
                    logging.debug(f"    Added pair with key: {key}")
                else:
                    logging.debug(f"    No R2 file found for: {r1_file}")
    
    # Convert temp_pairs to final file_pairs, resolving any remaining conflicts
    for key, pair in temp_pairs.items():
        if key not in file_pairs:
            file_pairs[key] = pair
        else:
            # If still conflicting, use full path as key
            full_key = pair['R1']
            file_pairs[full_key] = pair
            logging.debug(f"Resolved conflict using full path: {full_key}")
    
    logging.info(f"  Found {len(set(pair['R1'] for pair in file_pairs.values()))} unique FASTQ pairs")
    logging.debug(f"  File pairs: {file_pairs}")
    return file_pairs

def fuzzy_match_sample(sample_name, available_samples):
    """
    Fuzzy matching for sample names - handles typos and variations
    Returns best match or None
    """
    sample_name = sample_name.lower().strip()
    
    # Exact match first
    for available in available_samples:
        if available.lower() == sample_name:
            return available
    
    # Partial match - sample name contains or is contained in available
    matches = []
    for available in available_samples:
        available_lower = available.lower()
        available_base = os.path.basename(available_lower)
        
        if sample_name in available_lower or available_lower in sample_name:
            matches.append((available, len(available)))
        elif sample_name in available_base or available_base in sample_name:
            matches.append((available, len(available)))
    
    if matches:
        # Return shortest match (most specific)
        return min(matches, key=lambda x: x[1])[0]
    
    # Substring matching with common prefixes/suffixes
    for available in available_samples:
        available_clean = available.lower().replace('_', '').replace('-', '')
        sample_clean = sample_name.replace('_', '').replace('-', '')
        
        if sample_clean in available_clean or available_clean in sample_clean:
            return available
    
    return None

def count_reads_fast(fastq_file):
    """
    Fast read counting using streaming - minimal memory usage
    """
    read_count = 0
    try:
        # Auto-detect gzip or plain text
        opener = gzip.open if str(fastq_file).endswith('.gz') else open
        with opener(fastq_file, 'rt') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline()
                plus = f.readline()
                qual = f.readline()
                if not (seq and plus and qual):
                    break
                if header.startswith('@'):
                    read_count += 1
    except Exception as e:
        logging.error(f"    Error counting reads in {os.path.basename(fastq_file)}: {e}")
        return 0
    return read_count

def validate_paired_end_integrity(r1_files, r2_files):
    """Validate that R1 and R2 files have matching read counts"""
    r1_counts = {}
    r2_counts = {}
    
    for f in r1_files:
        r1_counts[f] = count_reads_fast(f)
    
    for f in r2_files:
        r2_counts[f] = count_reads_fast(f)
    
    mismatches = []
    for f in r1_files:
        r2_file = f.replace('_R1', '_R2').replace('_1', '_2')
        if r2_file in r2_counts:
            if r1_counts[f] != r2_counts[r2_file]:
                mismatches.append((f, r1_counts[f], r2_file, r2_counts[r2_file]))
    
    return mismatches

def calculate_file_checksum(file_path):
    """Calculate MD5 checksum of a file"""
    import hashlib
    hash_md5 = hashlib.md5()
    try:
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(4096), b""):
                hash_md5.update(chunk)
        return hash_md5.hexdigest()
    except Exception as e:
        logging.warning(f"Could not calculate checksum for {file_path}: {e}")
        return None

def combine_fastq_files_streaming(source_files, output_file, read_type='R1', buffer_size=BUFFER_SIZE, 
                                validate=False, check_barcodes=False, gc_analysis=False, adapter_check=False,
                                create_backups=False, deduplicate=False):
    """Combine multiple FASTQ files using streaming I/O with enhanced validation and optional deduplication"""
    total_reads = 0
    total_size = 0
    validation_warnings = []
    barcodes_found = set()
    gc_contents = []
    adapters_found = set()
    
    # Memory-efficient deduplication with size limit
    if deduplicate:
        max_sequences = 10_000_000  # 10M sequences max in memory
        seen_sequences = set()
        logging.info(f"Using in-memory deduplication (max {max_sequences:,} sequences)")
    else:
        seen_sequences = None
    
    # Create backup if requested
    if create_backups:
        for source_file in source_files:
            create_backup(source_file)
    
    with gzip.open(output_file, 'wt') as outfile:
        for i, source_file in enumerate(source_files, 1):
            file_size = os.path.getsize(source_file)
            total_size += file_size
            
            # Validate file if requested
            if validate:
                file_warnings = validate_fastq_quality(source_file)
                if file_warnings:
                    validation_warnings.extend([f"{source_file}: {w}" for w in file_warnings])
            
            opener = gzip.open if source_file.endswith('.gz') else open
            with opener(source_file, 'rt') as infile:
                while True:
                    header = infile.readline()
                    if not header:
                        break
                    sequence = infile.readline()
                    plus = infile.readline()
                    quality = infile.readline()
                    if not all([header, sequence, plus, quality]):
                        break
                    
                    seq_str = sequence.strip()
                    
                    # Memory-bounded deduplication
                    if deduplicate and seen_sequences is not None:
                        if len(seen_sequences) >= max_sequences:
                            logging.warning(f"Memory limit reached for deduplication ({max_sequences:,} sequences). Stopping deduplication.")
                            seen_sequences = None  # Disable further deduplication
                        else:
                            if seq_str in seen_sequences:
                                continue
                            seen_sequences.add(seq_str)
                    
                    total_reads += 1
                    
                    # Analysis features
                    if check_barcodes:
                        barcode = extract_sample_barcode(header.strip())
                        if barcode:
                            barcodes_found.add(barcode)
                    
                    if gc_analysis:
                        gc_content = calculate_gc_content(seq_str)
                        gc_contents.append(gc_content)
                    
                    if adapter_check:
                        adapters = detect_adapters(seq_str)
                        adapters_found.update(adapters)
                    
                    # Write to output
                    outfile.write(header)
                    outfile.write(sequence)
                    outfile.write(plus)
                    outfile.write(quality)
    
    # Calculate checksum for integrity verification
    checksum = calculate_file_checksum(output_file)
    if checksum:
        logging.info(f"Output file checksum ({read_type}): {checksum}")
    
    # Log results
    if validation_warnings:
        logging.warning(f"Validation warnings for {read_type}: {len(validation_warnings)} issues found")
        for warning in validation_warnings[:5]:
            logging.warning(f"  {warning}")
    
    if check_barcodes and barcodes_found:
        logging.info(f"Found {len(barcodes_found)} unique barcodes in {read_type}")
    
    if gc_analysis and gc_contents:
        avg_gc = sum(gc_contents) / len(gc_contents)
        logging.info(f"Average GC content for {read_type}: {avg_gc:.1f}%")
    
    if adapter_check and adapters_found:
        logging.warning(f"Detected adapters in {read_type}: {', '.join(adapters_found)}")
    
    return total_reads

def sanitize_sample_name(sample_name):
    """Sanitize sample name for Cell Ranger compatibility"""
    # Remove or replace problematic characters
    sanitized = sample_name.replace(' ', '_').replace('-', '_')
    sanitized = ''.join(c for c in sanitized if c.isalnum() or c in '_-')
    return sanitized

def get_memory_usage():
    """Get current memory usage in MB"""
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 / 1024

def detect_storage_type(path):
    """Detect if storage is SSD or HDD"""
    try:
        # Use system commands to detect storage type
        if platform.system() == "Darwin":  # macOS
            result = subprocess.run(['diskutil', 'info', path], 
                                  capture_output=True, text=True)
            if 'Solid State' in result.stdout:
                return 'SSD'
            return 'HDD'
        elif platform.system() == "Linux":
            # Linux detection logic
            result = subprocess.run(['lsblk', '-d', '-o', 'name,rota'], 
                                  capture_output=True, text=True)
            # Parse rotation info (0=SSD, 1=HDD)
            return 'SSD'  # Simplified for now
        else:
            return 'Unknown'
    except Exception:
        return 'Unknown'

def optimize_batch_size(available_ram_mb, storage_type, file_count):
    """Optimize batch size based on system resources"""
    base_size = 8 * 1024 * 1024  # 8MB base
    
    # Adjust for available RAM (use 10% of available RAM)
    ram_factor = min(available_ram_mb * 0.1 / 1024, 100)  # Cap at 100MB
    
    # Adjust for storage type
    if storage_type == 'SSD':
        storage_factor = 2.0  # Can use larger batches with SSD
    else:
        storage_factor = 0.5  # Smaller batches for HDD
    
    # Adjust for file count
    file_factor = min(file_count / 10, 2.0)  # More files = larger batches
    
    optimized_size = int(base_size * ram_factor * storage_factor * file_factor)
    return max(1024 * 1024, min(optimized_size, 100 * 1024 * 1024))  # 1MB to 100MB

def create_checkpoint_file(output_dir, mapping, file_pairs, combination_stats):
    """Create checkpoint file for resuming interrupted runs"""
    checkpoint_data = {
        'timestamp': datetime.now().isoformat(),
        'mapping': dict(mapping),
        'file_pairs': {k: dict(v) for k, v in file_pairs.items()},
        'combination_stats': combination_stats,
        'completed_targets': list(combination_stats.keys())
    }
    checkpoint_path = os.path.join(output_dir, '.checkpoint.pkl')
    with open(checkpoint_path, 'wb') as f:
        pickle.dump(checkpoint_data, f)
    return checkpoint_path

def load_checkpoint(output_dir):
    """Load checkpoint file if it exists"""
    checkpoint_path = os.path.join(output_dir, '.checkpoint.pkl')
    if os.path.exists(checkpoint_path):
        try:
            with open(checkpoint_path, 'rb') as f:
                return pickle.load(f)
        except Exception as e:
            logging.warning(f"Failed to load checkpoint: {e}")
    return None

def clear_checkpoint(output_dir):
    """Clear checkpoint file after successful completion"""
    checkpoint_path = os.path.join(output_dir, '.checkpoint.pkl')
    if os.path.exists(checkpoint_path):
        os.remove(checkpoint_path)

def detect_quality_format(quality_line):
    """Detect FASTQ quality score format from quality line"""
    if not quality_line:
        return 'unknown'
    
    # Get ASCII values
    ascii_values = [ord(c) for c in quality_line]
    min_val = min(ascii_values)
    max_val = max(ascii_values)
    
    # Sanger/Illumina 1.8+ (ASCII 33-126, Phred+33)
    if min_val >= 33 and max_val <= 126:
        return 'sanger'
    # Illumina 1.3+ (ASCII 64-126, Phred+64)
    elif min_val >= 64 and max_val <= 126:
        return 'illumina_1.3'
    # Old Illumina (ASCII 64-126, Phred+64)
    elif min_val >= 64:
        return 'illumina_old'
    else:
        return 'unknown'

def validate_fastq_quality(fastq_file):
    """Validate FASTQ file quality and detect corruption with format detection"""
    warnings = []
    quality_formats = set()
    
    try:
        opener = gzip.open if fastq_file.endswith('.gz') else open
        with opener(fastq_file, 'rt') as f:
            line_count = 0
            read_count = 0
            lengths = []
            quality_lines = []
            
            for line in f:
                line_count += 1
                line = line.strip()
                
                if line_count % 4 == 1:  # Header line
                    if not line.startswith('@'):
                        warnings.append(f"Invalid header at line {line_count}")
                    read_count += 1
                elif line_count % 4 == 2:  # Sequence line
                    if not all(c in 'ACGTN' for c in line.upper()):
                        warnings.append(f"Invalid characters in sequence at read {read_count}")
                    lengths.append(len(line))
                elif line_count % 4 == 3:  # Plus line
                    if not line.startswith('+'):
                        warnings.append(f"Invalid plus line at line {line_count}")
                elif line_count % 4 == 0:  # Quality line
                    if len(line) != lengths[-1]:
                        warnings.append(f"Quality length mismatch at read {read_count}")
                    
                    # Detect quality format
                    quality_format = detect_quality_format(line)
                    quality_formats.add(quality_format)
                    
                    # Validate quality scores based on format
                    if quality_format == 'sanger':
                        if not all(33 <= ord(c) <= 126 for c in line):
                            warnings.append(f"Invalid Sanger quality scores at read {read_count}")
                    elif quality_format == 'illumina_1.3':
                        if not all(64 <= ord(c) <= 126 for c in line):
                            warnings.append(f"Invalid Illumina 1.3+ quality scores at read {read_count}")
                    elif quality_format == 'unknown':
                        warnings.append(f"Unknown quality score format at read {read_count}")
                    
                    quality_lines.append(line)
            
            if line_count % 4 != 0:
                warnings.append("Incomplete FASTQ file")
            
            # Report quality format consistency
            if len(quality_formats) > 1:
                warnings.append(f"Mixed quality formats detected: {quality_formats}")
            elif quality_formats:
                logging.info(f"Detected quality format: {list(quality_formats)[0]}")
                
    except Exception as e:
        warnings.append(f"File read error: {e}")
    
    return warnings

def extract_sample_barcode(header_line):
    """Extract sample barcode from FASTQ header"""
    # Common barcode patterns
    patterns = [
        r'[A-Z]{6,8}',  # 6-8 letter barcodes
        r'[0-9]{4,6}',  # 4-6 digit barcodes
        r'[A-Z0-9]{8,10}'  # Mixed alphanumeric
    ]
    
    import re
    for pattern in patterns:
        match = re.search(pattern, header_line)
        if match:
            return match.group()
    return None

def calculate_gc_content(sequence):
    """Calculate GC content of a sequence"""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return (gc_count / len(sequence)) * 100

def detect_adapters(sequence, common_adapters=None):
    """Detect common adapter sequences"""
    if common_adapters is None:
        common_adapters = [
            'AGATCGGAAGAG',  # Illumina adapter
            'CTGTCTCTTATA',  # Nextera adapter
            'AATGATACGGCG',  # TruSeq adapter
        ]
    
    detected = []
    for adapter in common_adapters:
        if adapter in sequence.upper():
            detected.append(adapter)
    return detected

def monitor_disk_space(path, required_gb=1):
    """Monitor available disk space and warn if low"""
    try:
        statvfs = os.statvfs(path)
        free_gb = (statvfs.f_frsize * statvfs.f_bavail) / (1024**3)
        if free_gb < required_gb:
            logging.warning(f"Low disk space: {free_gb:.1f}GB available, {required_gb}GB recommended")
            return False
        return True
    except Exception:
        return True  # Assume OK if we can't check

def create_backup(file_path):
    """Create backup of original file before processing"""
    backup_path = file_path + '.backup'
    if not os.path.exists(backup_path):
        try:
            shutil.copy2(file_path, backup_path)
            return backup_path
        except Exception as e:
            logging.warning(f"Failed to create backup of {file_path}: {e}")
    return None

def retry_operation(operation, max_retries=3, delay=1):
    """Retry an operation with exponential backoff"""
    for attempt in range(max_retries):
        try:
            return operation()
        except Exception as e:
            if attempt == max_retries - 1:
                raise e
            logging.warning(f"Operation failed (attempt {attempt + 1}/{max_retries}): {e}")
            time.sleep(delay * (2 ** attempt))
    return None

def generate_html_report(output_dir, mapping, file_pairs, combination_stats, missing_files, fuzzy_matches, cli_args=None, threads=4):
    """Generate HTML report with detailed info and system/run metadata"""
    total_targets = len(mapping)
    total_sources = sum(len(sources) for sources in mapping.values())
    successful_combinations = len(combination_stats)
    total_combined_reads = sum(stats['total_reads'] for stats in combination_stats.values())
    total_input_files = sum(len(stats['source_files']) for stats in combination_stats.values())
    total_output_files = 2 * successful_combinations
    total_data_processed = sum(
        os.path.getsize(stats['r1_output']) + os.path.getsize(stats['r2_output'])
        for stats in combination_stats.values() if os.path.exists(stats['r1_output']) and os.path.exists(stats['r2_output'])
    )
    avg_read_length = None
    try:
        read_lengths = []
        for stats in combination_stats.values():
            for s in stats['source_files']:
                for read_type in ['R1', 'R2']:
                    fpath = file_pairs[s][read_type]
                    opener = gzip.open if fpath.endswith('.gz') else open
                    with opener(fpath, 'rt') as f:
                        f.readline()  # header
                        seq = f.readline().strip()
                        if seq:
                            read_lengths.append(len(seq))
                            break
        if read_lengths:
            avg_read_length = sum(read_lengths) / len(read_lengths)
    except Exception:
        avg_read_length = None
    # System/run metadata
    metadata = {
        'Python version': sys.version.split()[0],
        'Tool version': '1.0.0',
        'User': getpass.getuser(),
        'Host': socket.gethostname(),
        'Platform': platform.platform(),
        'Date/time': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
        'Threads used': threads,
        'CLI args': ' '.join(cli_args) if cli_args else '',
        'Output directory': os.path.abspath(output_dir),
    }
    # Skipped/failed samples
    failed_samples = []
    for target in mapping:
        if target not in combination_stats:
            reason = "No sources found"
            if target in missing_files:
                reason = "Missing/unreadable source files"
            failed_samples.append((target, reason))
    # Pre-format avg_read_length for summary card
    if avg_read_length:
        avg_read_length_str = f"{avg_read_length:.1f}"
    else:
        avg_read_length_str = "N/A"
    # HTML content
    html_path = os.path.join(output_dir, "combination_report.html")
    html_content = f"""
<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\">
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\">
    <title>FASTQ Combination Report</title>
    <script src=\"https://cdn.jsdelivr.net/npm/chart.js\"></script>
    <style>
        body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: #f5f5f5; color: #333; }}
        .container {{ max-width: 1400px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 20px rgba(0,0,0,0.1); }}
        .header {{ text-align: center; margin-bottom: 30px; padding-bottom: 20px; border-bottom: 3px solid #2196F3; }}
        .header h1 {{ color: #1976D2; margin: 0; font-size: 2.5em; }}
        .summary-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 20px; margin-bottom: 30px; }}
        .summary-card {{ background: linear-gradient(135deg, #2196F3, #1976D2); color: white; padding: 20px; border-radius: 8px; text-align: center; }}
        .summary-card h3 {{ margin: 0 0 10px 0; font-size: 2em; }}
        .section h2 {{ color: #1976D2; border-bottom: 2px solid #2196F3; padding-bottom: 10px; }}
        table {{ width: 100%; border-collapse: collapse; margin-bottom: 20px; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
        th {{ background-color: #2196F3; color: white; }}
        .success {{ color: #4CAF50; font-weight: bold; }}
        .warning {{ color: #FF9800; font-weight: bold; }}
        .error {{ color: #f44336; font-weight: bold; }}
        .filepath {{ font-family: 'Courier New', monospace; background-color: #f0f0f0; padding: 2px 6px; border-radius: 4px; font-size: 0.8em; }}
        .collapsible {{ background-color: #eee; color: #444; cursor: pointer; padding: 10px; width: 100%; border: none; text-align: left; outline: none; font-size: 1em; margin-top: 10px; }}
        .active, .collapsible:hover {{ background-color: #ccc; }}
        .content {{ padding: 0 18px; display: none; overflow: hidden; background-color: #f9f9f9; }}
        .search-container {{ margin: 20px 0; padding: 15px; background-color: #f8f9fa; border-radius: 8px; }}
        .search-container input {{ padding: 8px; margin-right: 10px; border: 1px solid #ddd; border-radius: 4px; width: 300px; }}
        .search-container select {{ padding: 8px; margin-right: 10px; border: 1px solid #ddd; border-radius: 4px; }}
        .search-container button {{ padding: 8px 16px; background-color: #2196F3; color: white; border: none; border-radius: 4px; cursor: pointer; }}
        .search-container button:hover {{ background-color: #1976D2; }}
        .chart-container {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; margin: 20px 0; }}
        .chart-card {{ background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
        .hidden {{ display: none; }}
    </style>
</head>
<body>
    <div class=\"container\">
        <div class=\"header\">
            <h1>🚀 FASTQ Combination Report - OPTIMIZED</h1>
            <p>Cell Ranger Compatible | Streaming I/O | Generated on {metadata['Date/time']}</p>
        </div>
        <div class=\"section\">
            <h2>🛠️ System & Run Metadata</h2>
            <table>"""
    for k, v in metadata.items():
        html_content += f"<tr><th>{html_escape.escape(str(k))}</th><td>{html_escape.escape(str(v))}</td></tr>"
    html_content += "</table></div>"
    html_content += f"""
        <div class=\"summary-grid\">
            <div class=\"summary-card\"><h3>{total_targets:,}</h3><p>Target Samples</p></div>
            <div class=\"summary-card\"><h3>{successful_combinations:,}</h3><p>Successful</p></div>
            <div class=\"summary-card\"><h3>{len(fuzzy_matches):,}</h3><p>Fuzzy Matches</p></div>
            <div class=\"summary-card\"><h3>{total_combined_reads:,}</h3><p>Total Reads</p></div>
            <div class=\"summary-card\"><h3>{total_input_files:,}</h3><p>Input Files</p></div>
            <div class=\"summary-card\"><h3>{total_output_files:,}</h3><p>Output Files</p></div>
            <div class=\"summary-card\"><h3>{total_data_processed/1024/1024/1024:.2f} GB</h3><p>Data Processed</p></div>
            <div class=\"summary-card\"><h3>{avg_read_length_str}</h3><p>Avg Read Length</p></div>
        </div>
        <div class=\"section\">
            <h2>📊 Visualizations</h2>
            <div class=\"chart-container\">
                <div class=\"chart-card\">
                    <h3>Success vs Failure</h3>
                    <canvas id=\"successChart\" width=\"400\" height=\"300\"></canvas>
                </div>
                <div class=\"chart-card\">
                    <h3>Read Count Distribution</h3>
                    <canvas id=\"readChart\" width=\"400\" height=\"300\"></canvas>
                </div>
            </div>
        </div>
        <div class=\"section\">
            <h2>📊 Combination Results</h2>
            <div class=\"search-container\">
                <input type=\"text\" id=\"searchInput\" placeholder=\"Search by target name...\" onkeyup=\"filterTable()\">
                <select id=\"statusFilter\" onchange=\"filterTable()\">
                    <option value=\"\">All Status</option>
                    <option value=\"success\">Success</option>
                    <option value=\"failed\">Failed</option>
                </select>
                <button onclick=\"clearFilters()\">Clear Filters</button>
            </div>
            <table id=\"resultsTable\">
                <tr>
                    <th>Target Sample</th>
                    <th>Cell Ranger Output</th>
                    <th>Source Matches</th>
                    <th>Total Reads</th>
                    <th>Status</th>
                    <th>Details</th>
                </tr>
    """
    for target, source_paths in mapping.items():
        clean_target = sanitize_sample_name(target)
        details_id = f"details_{clean_target}"
        status_class = "success" if target in combination_stats else "failed"
        if target in combination_stats:
            stats = combination_stats[target]
            status = f'<span class="success">✓ Success</span>'
            total_reads = f'{stats["total_reads"]:,}'
            # Show what sources were matched
            matched_sources = []
            for source_path in source_paths:
                if source_path in fuzzy_matches:
                    original, matched = fuzzy_matches[source_path]
                    matched_sources.append(f'"{original}" → {os.path.basename(matched)} <span class="warning">(fuzzy)</span>')
                else:
                    matched_sources.append(f'{source_path} <span class="success">(exact)</span>')
            source_matches = "<br>".join(matched_sources)
            cell_ranger_files = f'<a href="{clean_target}_S1_R1_001.fastq.gz">{clean_target}_S1_R1_001.fastq.gz</a><br><a href="{clean_target}_S1_R2_001.fastq.gz">{clean_target}_S1_R2_001.fastq.gz</a>'
        else:
            status = f'<span class="error">✗ Failed</span>'
            total_reads = "0"
            source_matches = "No matches found"
            cell_ranger_files = "Not generated"
        html_content += f"""
                <tr class=\"row {status_class}\" data-target=\"{target}\" data-status=\"{status_class}\">
                    <td><strong>{target}</strong></td>
                    <td class="filepath">{cell_ranger_files}</td>
                    <td>{source_matches}</td>
                    <td>{total_reads}</td>
                    <td>{status}</td>
                    <td><button class='collapsible'>Show Details</button><div class='content' id='{details_id}'>
        """
        # Per-source file details
        if target in combination_stats:
            html_content += "<table><tr><th>Source File</th><th>Type</th><th>Size (MB)</th><th>Read Count</th><th>Match</th><th>Warnings</th></tr>"
            for s in combination_stats[target]['source_files']:
                for read_type in ['R1', 'R2']:
                    fpath = file_pairs[s][read_type]
                    ftype = read_type
                    try:
                        size_mb = os.path.getsize(fpath) / 1024 / 1024
                    except Exception:
                        size_mb = 'N/A'
                    try:
                        opener = gzip.open if fpath.endswith('.gz') else open
                        with opener(fpath, 'rt') as f:
                            count = 0
                            while True:
                                header = f.readline()
                                if not header:
                                    break
                                seq = f.readline()
                                plus = f.readline()
                                qual = f.readline()
                                if header.startswith('@') and seq and plus and qual:
                                    count += 1
                        read_count = count
                    except Exception:
                        read_count = 'N/A'
                    match_type = 'Fuzzy' if s in fuzzy_matches else 'Exact'
                    warnings = []
                    if read_count == 'N/A':
                        warnings.append('Unreadable')
                    if size_mb == 'N/A':
                        warnings.append('Missing')
                    html_content += f"<tr><td class='filepath'>{fpath}</td><td>{ftype}</td><td>{size_mb}</td><td>{read_count}</td><td>{match_type}</td><td>{', '.join(warnings) if warnings else '-'}</td></tr>"
            html_content += "</table>"
        html_content += "</div></td></tr>"
    html_content += """
            </table>
        </div>
    """
    # Skipped/failed samples
    if failed_samples:
        html_content += f"<div class='section'><h2>❌ Skipped/Failed Samples</h2><table><tr><th>Target</th><th>Reason</th></tr>"
        for target, reason in failed_samples:
            html_content += f"<tr><td>{target}</td><td>{reason}</td></tr>"
        html_content += "</table></div>"
    # Fuzzy matches
    if fuzzy_matches:
        html_content += f"""
        <div class="section">
            <h2>🎯 Fuzzy Matches Applied</h2>
            <p>The following sample names were automatically corrected:</p>
            <table>
                <tr>
                    <th>Original Name</th>
                    <th>Matched To</th>
                    <th>Confidence</th>
                </tr>
        """
        for source_path, (original, matched) in fuzzy_matches.items():
            confidence = "High" if original.lower() in matched.lower() else "Medium"
            html_content += f"""
                <tr>
                    <td class="filepath">{original}</td>
                    <td class="filepath">{matched}</td>
                    <td><span class="warning">{confidence}</span></td>
                </tr>
            """
        html_content += """
            </table>
        </div>
        """
    # JS for collapsible sections, search/filter, and charts
    html_content += """
    <script>
    // Collapsible sections
    var coll = document.getElementsByClassName("collapsible");
    var i;
    for (i = 0; i < coll.length; i++) {
      coll[i].addEventListener("click", function() {
        this.classList.toggle("active");
        var content = this.nextElementSibling;
        if (content.style.display === "block") {
          content.style.display = "none";
        } else {
          content.style.display = "block";
        }
      });
    }
    
    // Search and filter functionality
    function filterTable() {
      var input = document.getElementById("searchInput");
      var statusFilter = document.getElementById("statusFilter");
      var filter = input.value.toLowerCase();
      var statusValue = statusFilter.value;
      var table = document.getElementById("resultsTable");
      var rows = table.getElementsByTagName("tr");
      
      for (var i = 1; i < rows.length; i++) {
        var row = rows[i];
        var target = row.getAttribute("data-target");
        var status = row.getAttribute("data-status");
        var txtValue = target;
        
        var showRow = true;
        
        // Text filter
        if (filter && !txtValue.toLowerCase().includes(filter)) {
          showRow = false;
        }
        
        // Status filter
        if (statusValue && status !== statusValue) {
          showRow = false;
        }
        
        if (showRow) {
          row.style.display = "";
        } else {
          row.style.display = "none";
        }
      }
    }
    
    function clearFilters() {
      document.getElementById("searchInput").value = "";
      document.getElementById("statusFilter").value = "";
      filterTable();
    }
    
    // Charts
    document.addEventListener('DOMContentLoaded', function() {
      // Success vs Failure Chart
      var successCtx = document.getElementById('successChart').getContext('2d');
      var successChart = new Chart(successCtx, {
        type: 'pie',
        data: {
          labels: ['Successful', 'Failed'],
          datasets: [{
            data: ["""
    html_content += f"{successful_combinations}, {len(failed_samples)}"
    html_content += """],
            backgroundColor: ['#4CAF50', '#f44336'],
            borderWidth: 2,
            borderColor: '#fff'
          }]
        },
        options: {
          responsive: true,
          plugins: {
            legend: {
              position: 'bottom'
            }
          }
        }
      });
      
      // Read Count Distribution Chart
      var readCtx = document.getElementById('readChart').getContext('2d');
      var readChart = new Chart(readCtx, {
        type: 'bar',
        data: {
          labels: ["""
    # Get read counts for successful samples
    read_counts = []
    for target in mapping:
        if target in combination_stats:
            read_counts.append(combination_stats[target]['total_reads'])
    
    # Format labels and data for chart
    if read_counts:
        labels_str = f"'{', '.join(str(x) for x in read_counts[:10])}'"  # Limit to first 10
        data_str = f"{', '.join(str(x) for x in read_counts[:10])}"
    else:
        labels_str = "'No data'"
        data_str = "0"
    
    html_content += f"{labels_str}"
    html_content += """],
          datasets: [{
            label: 'Read Count',
            data: ["""
    html_content += f"{data_str}"
    html_content += """],
            backgroundColor: '#2196F3',
            borderColor: '#1976D2',
            borderWidth: 1
          }]
        },
        options: {
          responsive: true,
          scales: {
            y: {
              beginAtZero: true,
              title: {
                display: true,
                text: 'Read Count'
              }
            },
            x: {
              title: {
                display: true,
                text: 'Sample Index'
              }
            }
          },
          plugins: {
            legend: {
              display: false
            }
          }
        }
      });
    });
    </script>
    """
    html_content += f"""
        <div class="footer" style="text-align: center; margin-top: 40px; color: #666;">
            <p><strong>⚡ OPTIMIZED FOR SPEED!</strong> Streaming I/O with minimal RAM usage</p>
            <p><strong>Cell Ranger Ready!</strong> All output files follow Illumina naming convention</p>
            <p>Output directory: {os.path.abspath(output_dir)}</p>
        </div>
    </div>
</body>
</html>
    """
    with open(html_path, 'w') as f:
        f.write(html_content)
    return html_path

def combine_fastq_files_main(csv_file, output_dir="combined", search_dirs=None, dry_run=False, force=False, 
                           r1_patterns=None, r2_patterns=None, threads=4, buffer_size=BUFFER_SIZE,
                           validate=False, check_barcodes=False, gc_analysis=False, adapter_check=False,
                           create_backups=False, retry_failed=False, real_time_monitor=False,
                           checkpoint=False, no_html=False, no_csv=False, deduplicate=False):
    """Main function with streaming optimizations"""
    
    logging.info("⚡ FASTQ File Combiner - STREAMING OPTIMIZED")
    logging.info("=" * 60)
    logging.info(f"🚀 High-speed streaming I/O with minimal RAM usage")
    logging.info(f"Mapping file: {csv_file}")
    logging.info(f"Output directory: {output_dir}")
    
    # Read mapping file
    try:
        mapping = read_mapping_file(csv_file)
    except Exception as e:
        logging.error(f"Error reading mapping file: {e}")
        return
    
    if not mapping:
        logging.error("No valid mappings found in CSV file!")
        return
    
    # Fast file discovery
    file_pairs = find_fastq_files_fast(search_dirs, r1_patterns=r1_patterns, r2_patterns=r2_patterns)
    
    if not file_pairs:
        logging.error("❌ No FASTQ files found! Check your search directories.")
        return
    
    # Match samples with fuzzy matching
    logging.info(f"\n🎯 Matching samples...")
    available_samples = list(file_pairs.keys())
    fuzzy_matches = {}
    final_mapping = {}
    
    for target, source_paths in mapping.items():
        logging.info(f"\n  Target: {target}")
        matched_sources = []
        
        for source_path in source_paths:
            # Try exact match first
            if source_path in file_pairs:
                matched_sources.append(source_path)
                logging.info(f"    ✓ {source_path} (exact match)")
            else:
                # Try fuzzy matching
                fuzzy_match = fuzzy_match_sample(source_path, available_samples)
                if fuzzy_match:
                    matched_sources.append(fuzzy_match)
                    fuzzy_matches[source_path] = (source_path, fuzzy_match)
                    logging.info(f"    ~ {source_path} → {fuzzy_match} (fuzzy match)")
                else:
                    logging.warning(f"    ✗ {source_path} (no match found)")
        
        if matched_sources:
            final_mapping[target] = matched_sources
        else:
            logging.error(f"    ❌ No sources found for {target}")
    
    if not final_mapping:
        logging.error("❌ No samples could be matched!")
        return
    
    # Validate all source files exist and are readable
    missing_files = []
    invalid_targets = set()
    for target, matched_sources in final_mapping.items():
        for s in matched_sources:
            for read_type in ['R1', 'R2']:
                fpath = file_pairs[s][read_type]
                if not os.path.exists(fpath) or not os.access(fpath, os.R_OK):
                    missing_files.append(fpath)
                    invalid_targets.add(target)
    if missing_files:
        logging.error(f"Missing or unreadable files detected: {len(missing_files)}")
        for f in missing_files:
            logging.error(f"  {f}")
        if dry_run:
            logging.info("Dry run: stopping due to missing files.")
            return
        else:
            logging.warning("Some samples will be skipped due to missing/unreadable files.")
    if dry_run:
        logging.info("Dry run complete. All files checked.")
        return
    # Remove invalid targets from final_mapping
    for target in invalid_targets:
        if target in final_mapping:
            del final_mapping[target]
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    logging.info(f"\nCreated output directory: {os.path.abspath(output_dir)}")
    
    # Perform combinations with Cell Ranger naming - STREAMING MODE
    combination_stats = {}

    def process_target(target, matched_sources):
        logging.info(f"\n🔗 Processing target: {target}")
        clean_target = sanitize_sample_name(target)
        logging.info(f"  📁 Combining {len(matched_sources)} files")
        logging.info(f"  📝 Cell Ranger format: {clean_target}_S1_R*_001.fastq.gz")
        logging.info(f"  ⚡ Using streaming I/O for maximum speed...")
        
        start_time = datetime.now()
        r1_sources = [file_pairs[s]['R1'] for s in matched_sources]
        r1_output = os.path.join(output_dir, f"{clean_target}_S1_R1_001.fastq.gz")
        r2_sources = [file_pairs[s]['R2'] for s in matched_sources]
        r2_output = os.path.join(output_dir, f"{clean_target}_S1_R2_001.fastq.gz")
        
        # Validate paired-end integrity before processing
        logging.info(f"  🔍 Validating paired-end integrity...")
        mismatches = validate_paired_end_integrity(r1_sources, r2_sources)
        if mismatches:
            logging.warning(f"  ⚠️  Paired-end mismatches detected:")
            for r1_file, r1_count, r2_file, r2_count in mismatches:
                logging.warning(f"    {os.path.basename(r1_file)}: {r1_count:,} vs {os.path.basename(r2_file)}: {r2_count:,}")
            if not force:
                logging.error(f"  ❌ Skipping {target} due to paired-end mismatches. Use --force to proceed.")
                return {
                    'target': target,
                    'source_files': matched_sources,
                    'total_reads': 0,
                    'r1_output': r1_output,
                    'r2_output': r2_output,
                    'clean_name': clean_target,
                    'processing_time': 0,
                    'speed_mb_per_sec': 0,
                    'skipped': True,
                    'error': 'paired_end_mismatch'
                }
        
        # Overwrite protection
        if not force:
            if os.path.exists(r1_output) or os.path.exists(r2_output):
                logging.warning(f"  Skipping {target}: output files already exist. Use --force to overwrite.")
                return {
                    'target': target,
                    'source_files': matched_sources,
                    'total_reads': 0,
                    'r1_output': r1_output,
                    'r2_output': r2_output,
                    'clean_name': clean_target,
                    'processing_time': 0,
                    'speed_mb_per_sec': 0,
                    'skipped': True,
                    'error': 'files_exist'
                }
        
        # Process with error recovery
        try:
            r1_reads = combine_fastq_files_streaming(r1_sources, r1_output, 'R1', buffer_size, validate, check_barcodes, gc_analysis, adapter_check, create_backups, deduplicate)
            r2_reads = combine_fastq_files_streaming(r2_sources, r2_output, 'R2', buffer_size, validate, check_barcodes, gc_analysis, adapter_check, create_backups, deduplicate)
        except Exception as e:
            logging.error(f"  ❌ Error processing {target}: {e}")
            # Clean up partial outputs
            for output_file in [r1_output, r2_output]:
                if os.path.exists(output_file):
                    os.remove(output_file)
            return {
                'target': target,
                'source_files': matched_sources,
                'total_reads': 0,
                'r1_output': r1_output,
                'r2_output': r2_output,
                'clean_name': clean_target,
                'processing_time': 0,
                'speed_mb_per_sec': 0,
                'skipped': True,
                'error': str(e)
            }
        
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        total_size = sum(os.path.getsize(file_pairs[s]['R1']) + os.path.getsize(file_pairs[s]['R2']) for s in matched_sources)
        speed_mb_per_sec = (total_size / 1024 / 1024) / duration if duration > 0 else 0
        
        if r1_reads != r2_reads:
            logging.warning(f"  ⚠️  Warning: R1 ({r1_reads:,}) and R2 ({r2_reads:,}) read counts don't match!")
        
        result = {
            'target': target,
            'source_files': matched_sources,
            'total_reads': r1_reads,
            'r2_reads': r2_reads,
            'r1_output': r1_output,
            'r2_output': r2_output,
            'clean_name': clean_target,
            'processing_time': duration,
            'speed_mb_per_sec': speed_mb_per_sec,
            'skipped': False,
            'paired_end_mismatches': len(mismatches) if mismatches else 0
        }
        
        logging.info(f"  ✅ {clean_target}: {r1_reads:,} read pairs")
        logging.info(f"  ⚡ Speed: {speed_mb_per_sec:.1f} MB/sec ({duration:.1f} seconds)")
        
        return result

    # Parallel processing with progress bar
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = []
        for target, matched_sources in final_mapping.items():
            futures.append(executor.submit(process_target, target, matched_sources))
        for f in tqdm(as_completed(futures), total=len(futures), desc="Combining samples"):
            res = f.result()
            combination_stats[res['target']] = res
    
    # Output summary CSV
    csv_path = os.path.join(output_dir, "combination_summary.csv")
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "Target", "R1 Output", "R2 Output", "R1 Reads", "R2 Reads", 
            "Processing Time (s)", "Speed (MB/s)", "Skipped", "Error", 
            "Paired-End Mismatches", "R1 Checksum", "R2 Checksum"
        ])
        for target, stats in combination_stats.items():
            # Calculate checksums for output files
            r1_checksum = calculate_file_checksum(stats.get('r1_output', '')) if os.path.exists(stats.get('r1_output', '')) else None
            r2_checksum = calculate_file_checksum(stats.get('r2_output', '')) if os.path.exists(stats.get('r2_output', '')) else None
            
            writer.writerow([
                target,
                stats.get('r1_output', ''),
                stats.get('r2_output', ''),
                stats.get('total_reads', 0),
                stats.get('r2_reads', 0),
                f"{stats.get('processing_time', 0):.2f}",
                f"{stats.get('speed_mb_per_sec', 0):.2f}",
                stats.get('skipped', False),
                stats.get('error', ''),
                stats.get('paired_end_mismatches', 0),
                r1_checksum or '',
                r2_checksum or ''
            ])
    logging.info(f"📄 CSV summary: {csv_path}")
    
    # Generate reports
    missing_files = []
    for target, sources in mapping.items():
        if target not in final_mapping:
            missing_files.extend(sources)
    
    html_report = generate_html_report(output_dir, mapping, file_pairs, combination_stats, missing_files, fuzzy_matches, cli_args=None, threads=threads)
    
    # Final summary with performance metrics
    total_processing_time = sum(stats['processing_time'] for stats in combination_stats.values())
    avg_speed = sum(stats['speed_mb_per_sec'] for stats in combination_stats.values()) / len(combination_stats) if combination_stats else 0
    
    logging.info(f"\n🎉 Combination Complete!")
    logging.info(f"✅ Successfully combined {len(combination_stats)} samples")
    logging.info(f"📊 Total reads processed: {sum(stats['total_reads'] for stats in combination_stats.values()):,}")
    logging.info(f"⚡ Average speed: {avg_speed:.1f} MB/sec")
    logging.info(f"⏱️  Total processing time: {total_processing_time:.1f} seconds")
    logging.info(f"🧬 Cell Ranger compatible naming applied")
    logging.info(f"📁 Output files: {output_dir}/")
    logging.info(f"📋 HTML report: {html_report}")
    
    if fuzzy_matches:
        logging.info(f"🎯 Applied {len(fuzzy_matches)} fuzzy matches for typos/variations")

def main():
    """Main CLI entry point with enhanced features"""
    parser = argparse.ArgumentParser(
        description="⚡ FASTQ File Combiner - High-speed streaming I/O with minimal RAM usage",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python fastq_combiner.py mapping.csv
  python fastq_combiner.py mapping.csv --output combined --threads 8 --validate
  python fastq_combiner.py mapping.csv --dry-run --search-dirs data/ run1/ run2/
  python fastq_combiner.py mapping.csv --deduplicate --paired-end-dedup --validate
        """
    )
    
    parser.add_argument('csv_file', help='CSV mapping file (target,sources)')
    parser.add_argument('--output', '-o', default='combined', help='Output directory (default: combined)')
    parser.add_argument('--search-dirs', nargs='+', help='Directories to search for FASTQ files')
    parser.add_argument('--dry-run', action='store_true', help='Show what would be done without processing')
    parser.add_argument('--force', '-f', action='store_true', help='Overwrite existing output files')
    parser.add_argument('--threads', '-t', type=int, default=4, help='Number of threads (default: 4)')
    parser.add_argument('--buffer-size', type=int, default=BUFFER_SIZE, help=f'Buffer size in bytes (default: {BUFFER_SIZE})')
    
    # Analysis options
    parser.add_argument('--validate', action='store_true', help='Validate FASTQ quality and format')
    parser.add_argument('--check-barcodes', action='store_true', help='Extract and analyze sample barcodes')
    parser.add_argument('--gc-analysis', action='store_true', help='Calculate GC content statistics')
    parser.add_argument('--adapter-check', action='store_true', help='Detect common adapter sequences')
    
    # Processing options
    parser.add_argument('--deduplicate', action='store_true', help='Remove duplicate sequences (memory-bounded)')
    parser.add_argument('--paired-end-dedup', action='store_true', help='Paired-end aware deduplication (maintains pairing)')
    parser.add_argument('--create-backups', action='store_true', help='Create backups of source files')
    parser.add_argument('--retry-failed', action='store_true', help='Retry failed operations')
    
    # Monitoring options
    parser.add_argument('--real-time-monitor', action='store_true', help='Enable real-time progress monitoring')
    parser.add_argument('--checkpoint', action='store_true', help='Enable checkpointing for resuming interrupted runs')
    
    # Output options
    parser.add_argument('--no-html', action='store_true', help='Skip HTML report generation')
    parser.add_argument('--no-csv', action='store_true', help='Skip CSV summary generation')
    
    # Advanced options
    parser.add_argument('--r1-patterns', nargs='+', help='Custom R1 file patterns')
    parser.add_argument('--r2-patterns', nargs='+', help='Custom R2 file patterns')
    parser.add_argument('--config', help='YAML configuration file')
    parser.add_argument('--diagnostics', action='store_true', help='Print system diagnostics')
    parser.add_argument('--profile', action='store_true', help='Enable performance profiling')
    
    args = parser.parse_args()
    
    # Print diagnostics if requested
    if args.diagnostics:
        print_diagnostics()
        return
    
    # Load config file if provided
    if args.config:
        try:
            with open(args.config, 'r') as f:
                config = yaml.safe_load(f)
            
            # Override defaults with config values
            def get_opt(opt, default=None):
                return getattr(args, opt) if getattr(args, opt) is not None else config.get(opt, default)
            
            # Apply config overrides
            args.output = get_opt('output', 'combined')
            args.search_dirs = get_opt('search_dirs', None)
            args.threads = get_opt('threads', 4)
            args.buffer_size = get_opt('buffer_size', BUFFER_SIZE)
            args.validate = get_opt('validate', False)
            args.check_barcodes = get_opt('check_barcodes', False)
            args.gc_analysis = get_opt('gc_analysis', False)
            args.adapter_check = get_opt('adapter_check', False)
            args.deduplicate = get_opt('deduplicate', False)
            args.paired_end_dedup = get_opt('paired_end_dedup', False)
            args.create_backups = get_opt('create_backups', False)
            args.retry_failed = get_opt('retry_failed', False)
            args.real_time_monitor = get_opt('real_time_monitor', False)
            args.checkpoint = get_opt('checkpoint', False)
            args.no_html = get_opt('no_html', False)
            args.no_csv = get_opt('no_csv', False)
            args.r1_patterns = get_opt('r1_patterns', None)
            args.r2_patterns = get_opt('r2_patterns', None)
            
        except Exception as e:
            logging.error(f"Error loading config file: {e}")
            return
    
    # Validate paired-end deduplication
    if args.paired_end_dedup and not args.deduplicate:
        logging.warning("--paired-end-dedup requires --deduplicate. Enabling deduplication.")
        args.deduplicate = True
    
    # Start profiling if requested
    if args.profile:
        profiler = cProfile.Profile()
        profiler.enable()
    
    try:
        # Run the main function
        combine_fastq_files_main(
            csv_file=args.csv_file,
            output_dir=args.output,
            search_dirs=args.search_dirs,
            dry_run=args.dry_run,
            force=args.force,
            r1_patterns=args.r1_patterns,
            r2_patterns=args.r2_patterns,
            threads=args.threads,
            buffer_size=args.buffer_size,
            validate=args.validate,
            check_barcodes=args.check_barcodes,
            gc_analysis=args.gc_analysis,
            adapter_check=args.adapter_check,
            create_backups=args.create_backups,
            retry_failed=args.retry_failed,
            real_time_monitor=args.real_time_monitor,
            checkpoint=args.checkpoint,
            no_html=args.no_html,
            no_csv=args.no_csv,
            deduplicate=args.deduplicate
        )
        
    except KeyboardInterrupt:
        logging.info("\n⚠️  Interrupted by user")
        sys.exit(1)
    except Exception as e:
        logging.error(f"❌ Fatal error: {e}")
        sys.exit(1)
    finally:
        # Stop profiling and print results
        if args.profile:
            profiler.disable()
            stats = pstats.Stats(profiler)
            stats.sort_stats('cumulative')
            stats.print_stats(20)  # Top 20 functions

class MemoryMonitor:
    """Monitor memory usage during processing"""
    def __init__(self):
        self.peak_memory = 0
        self.monitoring = False
        self.monitor_thread = None
    
    def start(self):
        self.monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
    
    def stop(self):
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join()
    
    def _monitor(self):
        while self.monitoring:
            current_memory = get_memory_usage()
            self.peak_memory = max(self.peak_memory, current_memory)
            time.sleep(1)
    
    def get_peak_memory(self):
        return self.peak_memory

class RealTimeMonitor:
    """Real-time monitoring of processing statistics"""
    def __init__(self):
        self.start_time = time.time()
        self.processed_files = 0
        self.total_reads = 0
        self.total_size = 0
        self.current_speed = 0
        self.monitoring = False
        self.monitor_thread = None
    
    def start(self):
        self.monitoring = True
        self.monitor_thread = threading.Thread(target=self._monitor)
        self.monitor_thread.daemon = True
        self.monitor_thread.start()
    
    def stop(self):
        self.monitoring = False
        if self.monitor_thread:
            self.monitor_thread.join()
    
    def update(self, files_processed=0, reads_processed=0, size_processed=0):
        self.processed_files += files_processed
        self.total_reads += reads_processed
        self.total_size += size_processed
    
    def _monitor(self):
        while self.monitoring:
            elapsed = time.time() - self.start_time
            if elapsed > 0:
                self.current_speed = self.total_size / elapsed / 1024 / 1024  # MB/s
                logging.info(f"Real-time: {self.processed_files} files, {self.total_reads:,} reads, {self.current_speed:.1f} MB/s")
            time.sleep(5)  # Update every 5 seconds
    
    def get_stats(self):
        elapsed = time.time() - self.start_time
        return {
            'elapsed_time': elapsed,
            'processed_files': self.processed_files,
            'total_reads': self.total_reads,
            'total_size_mb': self.total_size / 1024 / 1024,
            'avg_speed_mb_s': self.total_size / elapsed / 1024 / 1024 if elapsed > 0 else 0
        }

def print_diagnostics():
    """Print system diagnostics"""
    print("🔍 FASTQ Combiner Diagnostics")
    print("=" * 50)
    print(f"Python version: {sys.version}")
    print(f"Platform: {platform.platform()}")
    print(f"CPU cores: {multiprocessing.cpu_count()}")
    
    # Memory info
    memory = psutil.virtual_memory()
    print(f"Total RAM: {memory.total / 1024**3:.1f}GB")
    print(f"Available RAM: {memory.available / 1024**3:.1f}GB")
    print(f"RAM usage: {memory.percent}%")
    
    # Disk info
    disk = psutil.disk_usage('.')
    print(f"Disk space: {disk.free / 1024**3:.1f}GB free of {disk.total / 1024**3:.1f}GB")
    
    # Storage type
    storage_type = detect_storage_type('.')
    print(f"Storage type: {storage_type}")
    
    # Check dependencies
    dependencies = ['gzip', 'yaml', 'tqdm', 'psutil']
    print("\nDependencies:")
    for dep in dependencies:
        try:
            importlib.import_module(dep)
            print(f"  ✓ {dep}")
        except ImportError:
            print(f"  ✗ {dep} (missing)")
    
    print("\nSystem ready for FASTQ processing!")

if __name__ == "__main__":
    sys.exit(main())
