# Step 3: SAM to BAM Conversion, Sorting, and Indexing

## Objective
Convert SAM alignment files to the more efficient BAM format, sort them by genomic coordinates, and create index files for fast random access in downstream analysis.

## Why BAM Format is Essential

### SAM vs BAM:
- **SAM (Sequence Alignment/Map)**: Human-readable text format
- **BAM (Binary Alignment/Map)**: Compressed binary version of SAM

### Benefits of BAM:
- **File size**: 3-5x smaller than SAM files
- **Speed**: Much faster to read and process
- **Indexing**: Enables rapid access to specific genomic regions
- **Tool compatibility**: Required by most downstream ChIP-seq tools
- **Sorting**: Necessary for duplicate removal and peak calling

### Why Sorting is Critical:
- **Duplicate detection**: Picard/samtools need sorted reads to identify duplicates
- **Peak calling**: MACS2 requires sorted BAM files
- **Visualization**: Genome browsers need sorted, indexed BAM files
- **Memory efficiency**: Sorted files enable streaming algorithms

## Understanding Coordinate Sorting

Coordinate sorting arranges reads by:
1. **Chromosome**: chr1, chr2, ..., chr10
2. **Position**: Left-most mapping position within each chromosome
3. **Strand**: Forward strand before reverse strand at same position

This organization is crucial for ChIP-seq analysis because ARF27 binding sites need to be analyzed in genomic context.

## Step-by-Step Instructions

### 1. Verify SAMtools Installation

```bash
# Check if samtools is installed
samtools --version

# Install if needed
conda install -c bioconda samtools

# Verify version (recommended: ≥1.10)
samtools --version | head -1
```

### 2. Convert SAM to BAM

Convert all SAM files to BAM format:

```bash
# Navigate to your analysis directory
cd chipseq_arf27_analysis

# Convert ARF27 ChIP replicate 1
echo "Converting ARF27 ChIP replicate 1 to BAM..."
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam > 04_mapped_data/ARF27_ChIP_rep1_unsorted.bam

# Convert ARF27 ChIP replicate 2
echo "Converting ARF27 ChIP replicate 2 to BAM..."
samtools view -bS 04_mapped_data/ARF27_ChIP_rep2.sam > 04_mapped_data/ARF27_ChIP_rep2_unsorted.bam

# Convert Input control replicate 1
echo "Converting Input control replicate 1 to BAM..."
samtools view -bS 04_mapped_data/Input_control_rep1.sam > 04_mapped_data/Input_control_rep1_unsorted.bam

# Convert Input control replicate 2
echo "Converting Input control replicate 2 to BAM..."
samtools view -bS 04_mapped_data/Input_control_rep2.sam > 04_mapped_data/Input_control_rep2_unsorted.bam
```

**Command Explanation:**
- `samtools view`: Convert between SAM/BAM formats
- `-b`: Output in BAM format
- `-S`: Input is SAM format (auto-detected in newer versions)

### 3. Sort BAM Files by Coordinates

Sort BAM files for downstream analysis:

```bash
# Sort ARF27 ChIP replicate 1
echo "Sorting ARF27 ChIP replicate 1..."
samtools sort 04_mapped_data/ARF27_ChIP_rep1_unsorted.bam -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam

# Sort ARF27 ChIP replicate 2
echo "Sorting ARF27 ChIP replicate 2..."
samtools sort 04_mapped_data/ARF27_ChIP_rep2_unsorted.bam -o 04_mapped_data/ARF27_ChIP_rep2_sorted.bam

# Sort Input control replicate 1
echo "Sorting Input control replicate 1..."
samtools sort 04_mapped_data/Input_control_rep1_unsorted.bam -o 04_mapped_data/Input_control_rep1_sorted.bam

# Sort Input control replicate 2
echo "Sorting Input control replicate 2..."
samtools sort 04_mapped_data/Input_control_rep2_unsorted.bam -o 04_mapped_data/Input_control_rep2_sorted.bam
```

**Command Explanation:**
- `samtools sort`: Sort BAM file by coordinates
- `-o`: Specify output filename
- Default: Sorts by chromosome and position

### 4. Index Sorted BAM Files

Create index files for fast random access:

```bash
# Index ARF27 ChIP replicate 1
echo "Indexing ARF27 ChIP replicate 1..."
samtools index 04_mapped_data/ARF27_ChIP_rep1_sorted.bam

# Index ARF27 ChIP replicate 2
echo "Indexing ARF27 ChIP replicate 2..."
samtools index 04_mapped_data/ARF27_ChIP_rep2_sorted.bam

# Index Input control replicate 1
echo "Indexing Input control replicate 1..."
samtools index 04_mapped_data/Input_control_rep1_sorted.bam

# Index Input control replicate 2
echo "Indexing Input control replicate 2..."
samtools index 04_mapped_data/Input_control_rep2_sorted.bam
```

**Expected Output:**
Each BAM file gets a corresponding `.bai` index file:
- `ARF27_ChIP_rep1_sorted.bam.bai`
- `ARF27_ChIP_rep2_sorted.bam.bai`
- etc.

### 5. Alternative: One-Step Conversion, Sorting, and Indexing

For efficiency, you can combine these operations:

```bash
# One-step conversion and sorting (alternative approach)
# Use this instead of steps 2-4 if you prefer

# Process ARF27 ChIP replicate 1
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam | \
samtools sort -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam
samtools index 04_mapped_data/ARF27_ChIP_rep1_sorted.bam

# Process ARF27 ChIP replicate 2
samtools view -bS 04_mapped_data/ARF27_ChIP_rep2.sam | \
samtools sort -o 04_mapped_data/ARF27_ChIP_rep2_sorted.bam
samtools index 04_mapped_data/ARF27_ChIP_rep2_sorted.bam

# Process Input control replicate 1
samtools view -bS 04_mapped_data/Input_control_rep1.sam | \
samtools sort -o 04_mapped_data/Input_control_rep1_sorted.bam
samtools index 04_mapped_data/Input_control_rep1_sorted.bam

# Process Input control replicate 2
samtools view -bS 04_mapped_data/Input_control_rep2.sam | \
samtools sort -o 04_mapped_data/Input_control_rep2_sorted.bam
samtools index 04_mapped_data/Input_control_rep2_sorted.bam
```

### 6. Clean Up Intermediate Files

After confirming successful conversion:

```bash
# Check that all sorted BAM files and indices were created
ls -la 04_mapped_data/*sorted.bam*

# Verify BAM file integrity
samtools quickcheck 04_mapped_data/*sorted.bam

# If all looks good, remove intermediate files to save space
rm 04_mapped_data/*_unsorted.bam

# Optional: Remove SAM files if space is limited
# (Keep them if you have space for backup)
# rm 04_mapped_data/*.sam
```

## Quality Control and Verification

### 1. Verify File Integrity

```bash
# Check BAM file integrity
echo "Checking BAM file integrity..."
for bam in 04_mapped_data/*sorted.bam
do
    echo "Checking $bam"
    samtools quickcheck $bam
    if [ $? -eq 0 ]; then
        echo "✓ $bam is valid"
    else
        echo "✗ $bam is corrupted"
    fi
done
```

### 2. Compare File Sizes

```bash
# Compare file sizes
echo "File size comparison:"
ls -lh 04_mapped_data/ARF27_ChIP_rep1.sam 04_mapped_data/ARF27_ChIP_rep1_sorted.bam 2>/dev/null || echo "SAM files may have been removed"
ls -lh 04_mapped_data/*sorted.bam
```

**Expected size reduction**: BAM files should be 3-5x smaller than SAM files

### 3. Verify Read Counts

```bash
# Count reads in SAM vs BAM to ensure no loss
echo "Verifying read counts..."

# Count reads in original SAM file (if still available)
if [ -f "04_mapped_data/ARF27_ChIP_rep1.sam" ]; then
    sam_count=$(samtools view -c 04_mapped_data/ARF27_ChIP_rep1.sam)
    echo "SAM reads: $sam_count"
fi

# Count reads in sorted BAM file
bam_count=$(samtools view -c 04_mapped_data/ARF27_ChIP_rep1_sorted.bam)
echo "BAM reads: $bam_count"

# They should be identical
```

### 4. Test Random Access

```bash
# Test index functionality by extracting reads from a specific region
# Example: Extract reads from chromosome 1, positions 1-1000000
echo "Testing BAM index functionality..."
samtools view 04_mapped_data/ARF27_ChIP_rep1_sorted.bam chr1:1-1000000 | head -5

# This should return quickly if indexing worked properly
```

### 5. Generate Basic Statistics

```bash
# Generate comprehensive statistics for each BAM file
echo "Generating BAM statistics..."
for bam in 04_mapped_data/*sorted.bam
do
    sample=$(basename $bam .bam)
    echo "Statistics for $sample:"
    samtools flagstat $bam > logs/${sample}_flagstat.txt
    samtools idxstats $bam > logs/${sample}_idxstats.txt
    echo "✓ Statistics saved to logs/"
done
```

## Understanding the Output Files

### File Types Created:
1. **Sorted BAM files** (`*_sorted.bam`): Primary analysis files
2. **Index files** (`*.bam.bai`): Enable fast random access
3. **Statistics files** (`*_flagstat.txt`, `*_idxstats.txt`): QC metrics

### BAM File Properties:
- **Coordinate-sorted**: Reads arranged by genomic position
- **Compressed**: Binary format with efficient compression
- **Indexed**: Rapid access to any genomic region
- **Standard format**: Compatible with all ChIP-seq tools

## Memory and Performance Optimization

### For Large Files (>20GB):

```bash
# Use more memory for sorting (if available)
samtools sort -m 8G 04_mapped_data/ARF27_ChIP_rep1_unsorted.bam -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam

# Use multiple threads (if available)
samtools sort -@ 8 04_mapped_data/ARF27_ChIP_rep1_unsorted.bam -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam
```

### For Limited Disk Space:

```bash
# Stream processing to avoid temporary files
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam | \
samtools sort -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam -

# Remove SAM immediately after conversion
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam | \
samtools sort -o 04_mapped_data/ARF27_ChIP_rep1_sorted.bam - && \
rm 04_mapped_data/ARF27_ChIP_rep1.sam
```

## Common Issues and Troubleshooting

### Issue 1: "Truncated file" Error

**Symptoms:**
```
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated.
```

**Causes:**
- Incomplete file transfer
- Insufficient disk space during conversion
- Process was interrupted

**Solutions:**
```bash
# Check disk space
df -h

# Verify original SAM file integrity
samtools view -H 04_mapped_data/ARF27_ChIP_rep1.sam > /dev/null

# Re-run conversion with verbose output
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam > 04_mapped_data/ARF27_ChIP_rep1_unsorted.bam -v
```

### Issue 2: Sorting Takes Too Long

**Optimization strategies:**
```bash
# Monitor memory usage
htop

# Increase sort memory (default is 768M)
samtools sort -m 4G input.bam -o output_sorted.bam

# Use multiple threads
samtools sort -@ 4 input.bam -o output_sorted.bam

# Use temporary directory on fast storage
samtools sort -T /tmp/sort_temp input.bam -o output_sorted.bam
```

### Issue 3: Index Creation Fails

**Common causes:**
- BAM file is not sorted
- BAM file is corrupted
- Insufficient permissions

**Solutions:**
```bash
# Verify BAM is sorted
samtools view -H 04_mapped_data/ARF27_ChIP_rep1_sorted.bam | grep "SO:"

# Check file permissions
ls -la 04_mapped_data/

# Force re-sorting if needed
samtools sort 04_mapped_data/ARF27_ChIP_rep1_sorted.bam -o 04_mapped_data/ARF27_ChIP_rep1_resorted.bam
```

### Issue 4: High Memory Usage

**For systems with limited RAM:**
```bash
# Reduce sort memory
samtools sort -m 1G input.bam -o output_sorted.bam

# Process files one at a time instead of parallel
# Monitor with: htop or top
```

## File Organization Check

Verify your directory structure:

```bash
# Check final file organization
echo "Final file structure:"
ls -la 04_mapped_data/

# Expected files:
# - *_sorted.bam (main analysis files)
# - *.bam.bai (index files)
# - Optionally: *.sam (original alignment files)
# - logs/*_flagstat.txt (statistics)
# - logs/*_idxstats.txt (chromosome statistics)
```

## Expected File Sizes and Processing Times

### Typical File Sizes:
- **SAM files**: 20-80 GB per sample
- **BAM files**: 5-20 GB per sample (3-5x compression)
- **Index files**: 1-10 MB per BAM file

### Processing Times (8-core system):
- **SAM to BAM**: 10-30 minutes per sample
- **Sorting**: 15-45 minutes per sample
- **Indexing**: 1-5 minutes per sample

## Analysis Log Update

Document your conversion results:

```bash
# Update analysis log
cat >> logs/step3_conversion_log.txt << EOF

Step 3: SAM to BAM Conversion Results
Date: $(date)

Conversion Summary:
- All SAM files successfully converted to BAM
- All BAM files sorted by coordinates
- All BAM files indexed

File Sizes (BAM):
- ARF27_ChIP_rep1_sorted.bam: $(ls -lh 04_mapped_data/ARF27_ChIP_rep1_sorted.bam | awk '{print $5}')
- ARF27_ChIP_rep2_sorted.bam: $(ls -lh 04_mapped_data/ARF27_ChIP_rep2_sorted.bam | awk '{print $5}')
- Input_control_rep1_sorted.bam: $(ls -lh 04_mapped_data/Input_control_rep1_sorted.bam | awk '{print $5}')
- Input_control_rep2_sorted.bam: $(ls -lh 04_mapped_data/Input_control_rep2_sorted.bam | awk '{print $5}')

Quality Checks:
- BAM integrity: All files passed quickcheck
- Index functionality: Verified
- Read counts: Consistent with SAM files

Issues Encountered: [None/List any problems]

Ready for Step 4: PCR Duplicate Removal

EOF
```

## Key Takeaways

- BAM format is essential for efficient ChIP-seq analysis
- Coordinate sorting enables all downstream operations
- Index files provide rapid access to genomic regions
- Quality control at each step prevents downstream problems
- Proper file organization saves time and prevents errors

## Next Steps Checklist

Before proceeding to Step 4 (PCR Duplicate Removal):

- [ ] All BAM files successfully created and sorted
- [ ] All index files (.bai) present
- [ ] BAM integrity verified with quickcheck
- [ ] File sizes reasonable (3-5x smaller than SAM)
- [ ] Statistics files generated
- [ ] Documented any issues in analysis log

Ready for **Step 4: PCR Duplicate Removal**? We'll identify and remove PCR duplicates to improve the accuracy of ARF27 binding site detection.
