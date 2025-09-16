# Step 4: PCR Duplicate Removal for ChIP-seq Data

## Objective
Identify and remove PCR duplicates from ChIP-seq BAM files to improve signal-to-noise ratio and prevent false positive peak calls in ARF27 binding site detection.

## Why PCR Duplicate Removal is Critical for ChIP-seq

### What are PCR Duplicates?
PCR duplicates are multiple copies of the same DNA fragment that arise during library amplification. They appear as reads with:
- **Identical start positions** (5' coordinates)
- **Same fragment length** (for paired-end data)
- **Same chromosome and strand**

### Why Remove Duplicates in ChIP-seq?
1. **Artificial signal amplification**: Duplicates inflate read counts at binding sites
2. **False peak detection**: High duplicate regions can be mistaken for true binding
3. **Bias toward high-GC regions**: PCR tends to amplify GC-rich sequences preferentially
4. **Sequencing saturation**: Duplicates don't add new information
5. **Statistical assumptions**: Many tools assume independent observations

### ChIP-seq vs. Other Applications:
- **RNA-seq**: Duplicates often biological (high expression genes)
- **ChIP-seq**: Duplicates mostly technical (PCR artifacts)
- **Expected duplicate rate**: 15-40% for good ChIP-seq libraries

## Understanding Duplicate Detection

### Picard MarkDuplicates Algorithm:
1. **Group reads** by chromosome, strand, and 5' position
2. **For paired-end**: Also consider mate 5' position
3. **Within each group**: Keep the read with highest base quality sum
4. **Mark others**: Flag as duplicates but don't remove (option to remove)

### Special Considerations for ARF27 ChIP-seq:
- **Transcription factors** have focal binding sites
- **Higher duplicate rates** expected at true binding sites
- **Balance**: Remove technical duplicates while preserving biological signal

## Step-by-Step Instructions

### 1. Install Picard Tools

```bash
# Install Picard using conda (recommended)
conda install -c bioconda picard

# Or download from: https://broadinstitute.github.io/picard/
# Verify installation
picard MarkDuplicates --version
```

### 2. Create Output Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create directory for deduplicated BAM files
mkdir -p 05_processed_bams/deduplicated
```

### 3. Remove PCR Duplicates from All Samples

Process each sample individually:

#### ARF27 ChIP Replicate 1:

```bash
echo "Removing duplicates from ARF27 ChIP replicate 1..."
picard MarkDuplicates \
    INPUT=04_mapped_data/ARF27_ChIP_rep1_sorted.bam \
    OUTPUT=05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    METRICS_FILE=05_processed_bams/deduplicated/ARF27_ChIP_rep1_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    TMP_DIR=./tmp
```

#### ARF27 ChIP Replicate 2:

```bash
echo "Removing duplicates from ARF27 ChIP replicate 2..."
picard MarkDuplicates \
    INPUT=04_mapped_data/ARF27_ChIP_rep2_sorted.bam \
    OUTPUT=05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam \
    METRICS_FILE=05_processed_bams/deduplicated/ARF27_ChIP_rep2_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    TMP_DIR=./tmp
```

#### Input Control Replicate 1:

```bash
echo "Removing duplicates from Input control replicate 1..."
picard MarkDuplicates \
    INPUT=04_mapped_data/Input_control_rep1_sorted.bam \
    OUTPUT=05_processed_bams/deduplicated/Input_control_rep1_dedup.bam \
    METRICS_FILE=05_processed_bams/deduplicated/Input_control_rep1_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    TMP_DIR=./tmp
```

#### Input Control Replicate 2:

```bash
echo "Removing duplicates from Input control replicate 2..."
picard MarkDuplicates \
    INPUT=04_mapped_data/Input_control_rep2_sorted.bam \
    OUTPUT=05_processed_bams/deduplicated/Input_control_rep2_dedup.bam \
    METRICS_FILE=05_processed_bams/deduplicated/Input_control_rep2_dup_metrics.txt \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=LENIENT \
    CREATE_INDEX=true \
    TMP_DIR=./tmp
```

**Command Explanation:**
- `INPUT`: Sorted BAM file from Step 3
- `OUTPUT`: Deduplicated BAM file
- `METRICS_FILE`: Detailed duplication statistics
- `REMOVE_DUPLICATES=true`: Actually remove duplicates (vs. just marking)
- `ASSUME_SORTED=true`: Input is coordinate-sorted
- `VALIDATION_STRINGENCY=LENIENT`: Handle minor format issues
- `CREATE_INDEX=true`: Generate .bai index automatically
- `TMP_DIR`: Temporary directory for large file processing

### 4. Alternative: Mark Duplicates Without Removal

If you want to keep duplicates but mark them for exclusion:

```bash
# Mark duplicates without removing (alternative approach)
picard MarkDuplicates \
    INPUT=04_mapped_data/ARF27_ChIP_rep1_sorted.bam \
    OUTPUT=05_processed_bams/deduplicated/ARF27_ChIP_rep1_marked.bam \
    METRICS_FILE=05_processed_bams/deduplicated/ARF27_ChIP_rep1_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    CREATE_INDEX=true

# Later, exclude duplicates during analysis:
# samtools view -F 1024 input.bam  # Excludes duplicate flag
```

## Analyzing Duplicate Metrics

### 1. Examine Duplication Statistics

```bash
# View duplicate metrics for each sample
echo "=== ARF27 ChIP Rep 1 Duplicate Metrics ==="
cat 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dup_metrics.txt

echo "=== ARF27 ChIP Rep 2 Duplicate Metrics ==="
cat 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dup_metrics.txt

echo "=== Input Control Rep 1 Duplicate Metrics ==="
cat 05_processed_bams/deduplicated/Input_control_rep1_dup_metrics.txt

echo "=== Input Control Rep 2 Duplicate Metrics ==="
cat 05_processed_bams/deduplicated/Input_control_rep2_dup_metrics.txt
```

### 2. Extract Key Metrics

```bash
# Create summary of duplicate rates
echo "Sample\tTotal_Reads\tDuplicate_Reads\tPercent_Duplication" > logs/duplicate_summary.txt

for metrics in 05_processed_bams/deduplicated/*_dup_metrics.txt
do
    sample=$(basename $metrics _dup_metrics.txt)
    
    # Extract key metrics from Picard output
    total_reads=$(grep -A 1 "## METRICS CLASS" $metrics | tail -1 | cut -f3)
    duplicate_reads=$(grep -A 1 "## METRICS CLASS" $metrics | tail -1 | cut -f7)
    percent_dup=$(grep -A 1 "## METRICS CLASS" $metrics | tail -1 | cut -f9)
    
    echo -e "${sample}\t${total_reads}\t${duplicate_reads}\t${percent_dup}" >> logs/duplicate_summary.txt
done

# Display summary
column -t logs/duplicate_summary.txt
```

### 3. Understanding Metrics Output

Key metrics to examine:

#### LIBRARY (Main metrics line):
- **READ_PAIRS_EXAMINED**: Total read pairs processed
- **READ_PAIR_DUPLICATES**: Number of duplicate pairs
- **PERCENT_DUPLICATION**: Percentage of reads that are duplicates
- **ESTIMATED_LIBRARY_SIZE**: Estimated complexity of original library

#### Quality Indicators:
- **Good duplication rate**: 15-40% for ChIP-seq
- **Warning signs**: >60% duplication (over-amplified library)
- **Very low duplication**: <5% (possible under-amplification or quality issues)

## Quality Control and Validation

### 1. Compare Read Counts Before/After

```bash
echo "Read count comparison (before/after deduplication):"
echo "Sample\tOriginal_Reads\tAfter_Dedup\tRemoved\tPercent_Removed"

for sample in ARF27_ChIP_rep1 ARF27_ChIP_rep2 Input_control_rep1 Input_control_rep2
do
    original=$(samtools view -c 04_mapped_data/${sample}_sorted.bam)
    dedup=$(samtools view -c 05_processed_bams/deduplicated/${sample}_dedup.bam)
    removed=$((original - dedup))
    percent_removed=$(echo "scale=2; $removed * 100 / $original" | bc)
    
    echo -e "${sample}\t${original}\t${dedup}\t${removed}\t${percent_removed}%"
done
```

### 2. Verify File Integrity

```bash
# Check BAM file integrity
echo "Verifying deduplicated BAM file integrity..."
for bam in 05_processed_bams/deduplicated/*.bam
do
    echo "Checking $(basename $bam)"
    samtools quickcheck $bam
    if [ $? -eq 0 ]; then
        echo "✓ $(basename $bam) is valid"
    else
        echo "✗ $(basename $bam) is corrupted"
    fi
done
```

### 3. Generate Post-Deduplication Statistics

```bash
# Generate flagstat for deduplicated files
for bam in 05_processed_bams/deduplicated/*_dedup.bam
do
    sample=$(basename $bam _dedup.bam)
    echo "Generating statistics for $sample"
    samtools flagstat $bam > logs/${sample}_dedup_flagstat.txt
done
```

## Interpreting Duplication Results

### Expected Duplication Rates:

#### Good Quality Libraries:
- **ChIP samples**: 20-40% duplication
- **Input controls**: 15-30% duplication
- **Difference**: ChIP may have higher rates due to enrichment

#### Warning Signs:
- **Very high duplication (>60%)**:
  - Over-amplified library
  - Low starting material
  - Poor library complexity
  
- **Very low duplication (<5%)**:
  - Under-amplified library
  - Possible contamination
  - Technical issues

#### ARF27-Specific Considerations:
- **Transcription factors** typically have moderate duplication
- **Focal binding** may lead to localized high duplication
- **Compare replicates**: Should be similar within ~10%

### Sample Quality Assessment:

```bash
# Quick quality assessment
echo "Duplicate Rate Quality Assessment:"
echo "Sample\tDuplication_Rate\tQuality_Assessment"

for metrics in 05_processed_bams/deduplicated/*_dup_metrics.txt
do
    sample=$(basename $metrics _dup_metrics.txt)
    percent_dup=$(grep -A 1 "## METRICS CLASS" $metrics | tail -1 | cut -f9)
    
    # Convert to number for comparison
    dup_num=$(echo $percent_dup | sed 's/%//')
    
    if (( $(echo "$dup_num < 5" | bc -l) )); then
        quality="Low duplication - check library"
    elif (( $(echo "$dup_num <= 40" | bc -l) )); then
        quality="Good"
    elif (( $(echo "$dup_num <= 60" | bc -l) )); then
        quality="Moderate - acceptable"
    else
        quality="High - investigate library quality"
    fi
    
    echo -e "${sample}\t${percent_dup}\t${quality}"
done
```

## Common Issues and Troubleshooting

### Issue 1: Very High Duplication Rates (>60%)

**Possible Causes:**
- Over-amplification during PCR
- Low input DNA concentration
- Poor library complexity
- Sequencing saturation

**Diagnosis:**
```bash
# Check estimated library size from metrics
grep "ESTIMATED_LIBRARY_SIZE" 05_processed_bams/deduplicated/*_dup_metrics.txt

# Compare duplication rates between replicates
# Large differences may indicate technical problems
```

**Solutions:**
- If consistent across replicates: Proceed with analysis but note in methods
- If only one replicate affected: Consider excluding problematic replicate
- For future experiments: Optimize PCR conditions, increase input DNA

### Issue 2: Picard Out of Memory Error

**Error message:**
```
java.lang.OutOfMemoryError: Java heap space
```

**Solutions:**
```bash
# Increase Java memory allocation
export _JAVA_OPTIONS="-Xmx16g"

# Or specify memory in Picard command
picard -Xmx16g MarkDuplicates \
    INPUT=input.bam \
    OUTPUT=output.bam \
    METRICS_FILE=metrics.txt \
    REMOVE_DUPLICATES=true

# For very large files, use alternative tools
samtools rmdup input_sorted.bam output_rmdup.bam
```

### Issue 3: Different Duplication Rates Between Replicates

**Expected difference**: <10% between biological replicates
**Large differences (>20%)**:
- Check library preparation consistency
- Examine individual metrics files
- Consider technical factors

```bash
# Compare replicates
echo "Replicate comparison:"
grep -A 1 "## METRICS CLASS" 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dup_metrics.txt | tail -1 | cut -f9
grep -A 1 "## METRICS CLASS" 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dup_metrics.txt | tail -1 | cut -f9
```

### Issue 4: Temporary Directory Full

**Error**: No space left on device in tmp directory

**Solutions:**
```bash
# Use different temporary directory
mkdir -p /path/to/large/disk/tmp

# Specify in Picard command
picard MarkDuplicates \
    ... \
    TMP_DIR=/path/to/large/disk/tmp

# Clean up temporary files
rm -rf ./tmp/*
```

### Issue 5: Low Library Complexity

**Indicators:**
- High duplication rate
- Low estimated library size
- Plateau in duplication rate

**Assessment:**
```bash
# Check estimated library size
echo "Estimated library complexity:"
for metrics in 05_processed_bams/deduplicated/*_dup_metrics.txt
do
    sample=$(basename $metrics _dup_metrics.txt)
    lib_size=$(grep -A 1 "## METRICS CLASS" $metrics | tail -1 | cut -f10)
    echo "$sample: $lib_size"
done
```

## File Organization and Cleanup

### 1. Organize Output Files

```bash
# Check final file structure
ls -la 05_processed_bams/deduplicated/

# Expected files:
# - *_dedup.bam (deduplicated reads)
# - *_dedup.bai (index files)
# - *_dup_metrics.txt (duplication statistics)
```

### 2. Archive or Remove Intermediate Files

```bash
# Optional: Remove original sorted BAM files if space is limited
# (Only after confirming deduplication was successful)
# mv 04_mapped_data/*_sorted.bam archive/
# mv 04_mapped_data/*_sorted.bam.bai archive/
```

## Expected File Sizes and Processing Times

### File Size Changes:
- **Before deduplication**: 5-20 GB per BAM file
- **After deduplication**: 3-15 GB per BAM file (15-40% reduction)
- **Metrics files**: <1 MB each

### Processing Times (8-core system):
- **Small files (<5GB)**: 5-15 minutes
- **Medium files (5-15GB)**: 15-45 minutes  
- **Large files (>15GB)**: 45-90 minutes

## Analysis Log Update

```bash
# Update analysis log with duplication results
cat >> logs/step4_deduplication_log.txt << EOF

Step 4: PCR Duplicate Removal Results
Date: $(date)

Duplication Summary:
$(cat logs/duplicate_summary.txt)

Quality Assessment:
- ARF27_ChIP_rep1: $(grep ARF27_ChIP_rep1 logs/duplicate_summary.txt | cut -f4) duplication
- ARF27_ChIP_rep2: $(grep ARF27_ChIP_rep2 logs/duplicate_summary.txt | cut -f4) duplication  
- Input_control_rep1: $(grep Input_control_rep1 logs/duplicate_summary.txt | cut -f4) duplication
- Input_control_rep2: $(grep Input_control_rep2 logs/duplicate_summary.txt | cut -f4) duplication

Quality Checks:
- All BAM files passed integrity check
- Duplication rates within expected range: [Yes/No]
- Replicates show consistent duplication: [Yes/No]

Issues Encountered: [None/List problems and solutions]

Ready for Step 5: Data Normalization

EOF
```

## Key Takeaways

- PCR duplicate removal is essential for accurate ChIP-seq analysis
- Expected duplication rates: 15-40% for good ChIP-seq libraries
- High-quality libraries show consistent duplication between replicates
- ARF27 ChIP-seq may have focal high duplication at true binding sites
- Quality control at this step prevents downstream analysis artifacts

## Next Steps Checklist

Before proceeding to Step 5 (Data Normalization):

- [ ] All samples processed through duplicate removal
- [ ] Duplication rates documented and within expected range
- [ ] BAM file integrity verified
- [ ] Index files (.bai) generated
- [ ] Metrics files saved and analyzed
- [ ] Any quality issues documented and addressed

Ready for **Step 5: Data Normalization**? We'll normalize the deduplicated data using standard library size normalization and potentially spike-in adjustments to enable quantitative comparisons between samples and conditions.
