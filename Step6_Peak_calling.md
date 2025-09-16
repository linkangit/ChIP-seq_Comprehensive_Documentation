# Step 6: Peak Calling with MACS2 for ARF27 ChIP-seq

## Objective
Identify ARF27 binding sites in the maize genome by calling peaks from ChIP-seq data using MACS2, optimized for transcription factor ChIP-seq with appropriate statistical thresholds and control comparisons.

## Why Peak Calling is the Core of ChIP-seq Analysis

### What are ChIP-seq Peaks?
**Peaks** are genomic regions with significantly higher read density in ChIP samples compared to input controls, indicating protein binding sites.

### Peak Calling Challenges for Transcription Factors:
- **Sharp, narrow peaks**: TFs bind specific DNA motifs (50-200 bp)
- **Lower signal-to-noise**: Fewer binding sites than histone marks
- **Variable peak heights**: Different binding affinities
- **Background noise**: Non-specific binding and technical artifacts

### ARF27-Specific Considerations:
- **Auxin Response Factor**: Binds to Auxin Response Elements (AREs)
- **Expected motif**: TGTCTC or similar sequences
- **Binding pattern**: Sharp peaks at promoters and regulatory regions
- **Peak number**: Hundreds to thousands of binding sites expected

## Understanding MACS2 Algorithm

### MACS2 Peak Calling Process:
1. **Fragment size estimation**: Calculate average DNA fragment length
2. **Peak scanning**: Identify enriched regions using sliding windows
3. **Background modeling**: Build local λ (lambda) background model
4. **Statistical testing**: Poisson test for significance
5. **FDR correction**: Control false discovery rate
6. **Peak refinement**: Determine precise peak boundaries

### Key MACS2 Parameters for TF ChIP-seq:
- **--nomodel**: Skip shifting model for sharp TF peaks
- **--extsize**: Extend reads to fragment size
- **--shift**: Adjust for TF binding offset
- **--qvalue**: FDR threshold (typically 0.01-0.05)
- **--call-summits**: Find precise binding sites within peaks

## Step-by-Step Instructions

### 1. Install MACS2

```bash
# Install MACS2 using conda
conda install -c bioconda macs2

# Verify installation
macs2 --version

# Should show version 2.2.7 or newer
```

### 2. Prepare Peak Calling Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create peak calling directory
mkdir -p 06_peaks
mkdir -p 06_peaks/individual_replicates
mkdir -p 06_peaks/merged_replicates
mkdir -p 06_peaks/qc_plots
```

### 3. Call Peaks on Individual Replicates

First, call peaks on each replicate individually to assess consistency:

#### ARF27 ChIP Replicate 1:

```bash
echo "Calling peaks for ARF27 ChIP replicate 1..."
macs2 callpeak \
    --treatment 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    --control 05_processed_bams/deduplicated/Input_control_rep1_dedup.bam \
    --format BAM \
    --gsize 2.1e9 \
    --name ARF27_ChIP_rep1 \
    --outdir 06_peaks/individual_replicates \
    --qvalue 0.05 \
    --nomodel \
    --extsize 150 \
    --shift -75 \
    --call-summits \
    --keep-dup all \
    --verbose 3 \
    2> logs/ARF27_ChIP_rep1_peak_calling.log
```

#### ARF27 ChIP Replicate 2:

```bash
echo "Calling peaks for ARF27 ChIP replicate 2..."
macs2 callpeak \
    --treatment 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam \
    --control 05_processed_bams/deduplicated/Input_control_rep2_dedup.bam \
    --format BAM \
    --gsize 2.1e9 \
    --name ARF27_ChIP_rep2 \
    --outdir 06_peaks/individual_replicates \
    --qvalue 0.05 \
    --nomodel \
    --extsize 150 \
    --shift -75 \
    --call-summits \
    --keep-dup all \
    --verbose 3 \
    2> logs/ARF27_ChIP_rep2_peak_calling.log
```

**Command Explanation:**
- `--treatment`: ChIP-seq BAM file
- `--control`: Input control BAM file
- `--format BAM`: Input file format
- `--gsize 2.1e9`: Maize effective genome size (~2.1 Gb)
- `--name`: Output file prefix
- `--outdir`: Output directory
- `--qvalue 0.05`: FDR threshold (5%)
- `--nomodel`: Don't build shifting model (good for TFs)
- `--extsize 150`: Extend reads to 150bp (typical fragment size)
- `--shift -75`: Shift reads by half extension size
- `--call-summits`: Call precise peak summits
- `--keep-dup all`: Keep all duplicates (we already removed PCR dups)
- `--verbose 3`: Detailed logging

### 4. Call Peaks on Merged Replicates

For final peak set, merge replicates and call peaks:

```bash
# Merge ChIP replicates
echo "Merging ChIP replicates..."
samtools merge 06_peaks/merged_replicates/ARF27_ChIP_merged.bam \
    05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam

# Merge Input replicates
echo "Merging Input controls..."
samtools merge 06_peaks/merged_replicates/Input_control_merged.bam \
    05_processed_bams/deduplicated/Input_control_rep1_dedup.bam \
    05_processed_bams/deduplicated/Input_control_rep2_dedup.bam

# Index merged BAM files
samtools index 06_peaks/merged_replicates/ARF27_ChIP_merged.bam
samtools index 06_peaks/merged_replicates/Input_control_merged.bam

# Call peaks on merged data
echo "Calling peaks on merged replicates..."
macs2 callpeak \
    --treatment 06_peaks/merged_replicates/ARF27_ChIP_merged.bam \
    --control 06_peaks/merged_replicates/Input_control_merged.bam \
    --format BAM \
    --gsize 2.1e9 \
    --name ARF27_ChIP_merged \
    --outdir 06_peaks/merged_replicates \
    --qvalue 0.05 \
    --nomodel \
    --extsize 150 \
    --shift -75 \
    --call-summits \
    --keep-dup all \
    --verbose 3 \
    2> logs/ARF27_ChIP_merged_peak_calling.log
```

### 5. Alternative: Conservative Peak Calling

For high-confidence peaks, use more stringent parameters:

```bash
# Call peaks with stricter FDR threshold
echo "Calling high-confidence peaks..."
macs2 callpeak \
    --treatment 06_peaks/merged_replicates/ARF27_ChIP_merged.bam \
    --control 06_peaks/merged_replicates/Input_control_merged.bam \
    --format BAM \
    --gsize 2.1e9 \
    --name ARF27_ChIP_merged_strict \
    --outdir 06_peaks/merged_replicates \
    --qvalue 0.01 \
    --pvalue 1e-5 \
    --nomodel \
    --extsize 150 \
    --shift -75 \
    --call-summits \
    --keep-dup all \
    --verbose 3 \
    2> logs/ARF27_ChIP_merged_strict_peak_calling.log
```

### 6. Examine Peak Calling Output

MACS2 generates several output files:

```bash
# Check output files for merged replicate analysis
ls -la 06_peaks/merged_replicates/

# Expected files:
# ARF27_ChIP_merged_peaks.narrowPeak  - Main peak file (BED format)
# ARF27_ChIP_merged_peaks.xls         - Detailed peak information
# ARF27_ChIP_merged_summits.bed       - Peak summits (precise binding sites)
# ARF27_ChIP_merged_model.r           - R script for visualization
# ARF27_ChIP_merged_treat_pileup.bdg  - Treatment pileup track
# ARF27_ChIP_merged_control_lambda.bdg - Control lambda track
```

### 7. Basic Peak Statistics

```bash
# Count peaks called in each analysis
echo "Peak calling summary:"
echo "Sample\tTotal_Peaks\tStrong_Peaks"

for peaks in 06_peaks/individual_replicates/*_peaks.narrowPeak 06_peaks/merged_replicates/*_peaks.narrowPeak
do
    sample=$(basename $peaks _peaks.narrowPeak)
    total_peaks=$(wc -l < $peaks)
    # Strong peaks (score > 100, roughly q-value < 0.01)
    strong_peaks=$(awk '$5 > 100' $peaks | wc -l)
    echo -e "${sample}\t${total_peaks}\t${strong_peaks}"
done
```

## Quality Control and Peak Assessment

### 1. Analyze Peak Characteristics

```bash
# Examine peak lengths
echo "Analyzing peak lengths for merged dataset..."
awk '{print $3-$2}' 06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak > 06_peaks/qc_plots/peak_lengths.txt

# Calculate length statistics
R --vanilla << 'EOF'
lengths <- read.table("06_peaks/qc_plots/peak_lengths.txt", header=FALSE)$V1

cat("Peak length statistics:\n")
cat("Mean length:", round(mean(lengths), 1), "bp\n")
cat("Median length:", round(median(lengths), 1), "bp\n")
cat("Standard deviation:", round(sd(lengths), 1), "bp\n")
cat("Range:", min(lengths), "-", max(lengths), "bp\n")

# Create histogram
png("06_peaks/qc_plots/peak_length_distribution.png", width=800, height=600)
hist(lengths, breaks=50, 
     main="ARF27 Peak Length Distribution", 
     xlab="Peak Length (bp)", 
     ylab="Frequency",
     col="lightblue")
abline(v=median(lengths), col="red", lwd=2, lty=2)
legend("topright", paste("Median =", round(median(lengths)), "bp"), 
       col="red", lty=2, lwd=2)
dev.off()

# Expected: Most TF peaks 100-500bp, median ~200-300bp
EOF
```

### 2. Examine Peak Score Distribution

```bash
# Analyze peak scores (MACS2 score = -10*log10(qvalue))
echo "Analyzing peak scores..."
awk '{print $5}' 06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak > 06_peaks/qc_plots/peak_scores.txt

R --vanilla << 'EOF'
scores <- read.table("06_peaks/qc_plots/peak_scores.txt", header=FALSE)$V1

cat("Peak score statistics:\n")
cat("Mean score:", round(mean(scores), 1), "\n")
cat("Median score:", round(median(scores), 1), "\n")
cat("Number of high-confidence peaks (score > 100):", sum(scores > 100), "\n")

# Create score distribution plot
png("06_peaks/qc_plots/peak_score_distribution.png", width=800, height=600)
hist(scores, breaks=50, 
     main="ARF27 Peak Score Distribution", 
     xlab="MACS2 Score (-10*log10(qvalue))", 
     ylab="Frequency",
     col="lightgreen")
abline(v=100, col="red", lwd=2, lty=2)
text(100, max(hist(scores, plot=FALSE)$counts)/2, "q=0.01", pos=4, col="red")
dev.off()
EOF
```

### 3. Chromosome Distribution of Peaks

```bash
# Count peaks per chromosome
echo "Peaks per chromosome:"
echo "Chromosome\tPeak_Count\tPeaks_per_Mb"

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
do
    peak_count=$(awk -v chr="$chr" '$1==chr' 06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak | wc -l)
    
    # Get chromosome length (approximate)
    case $chr in
        chr1) chr_length=307 ;;
        chr2) chr_length=245 ;;
        chr3) chr_length=236 ;;
        chr4) chr_length=247 ;;
        chr5) chr_length=218 ;;
        chr6) chr_length=174 ;;
        chr7) chr_length=182 ;;
        chr8) chr_length=176 ;;
        chr9) chr_length=162 ;;
        chr10) chr_length=151 ;;
    esac
    
    peaks_per_mb=$(echo "scale=2; $peak_count / $chr_length" | bc)
    echo -e "${chr}\t${peak_count}\t${peaks_per_mb}"
done
```

### 4. Compare Replicates

```bash
# Compare peak overlap between replicates
echo "Comparing peak overlap between replicates..."

# Use bedtools to find overlapping peaks
bedtools intersect \
    -a 06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak \
    -b 06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak \
    -wa > 06_peaks/qc_plots/rep1_overlapping_rep2.bed

bedtools intersect \
    -a 06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak \
    -b 06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak \
    -wa > 06_peaks/qc_plots/rep2_overlapping_rep1.bed

# Calculate overlap statistics
rep1_total=$(wc -l < 06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak)
rep2_total=$(wc -l < 06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak)
rep1_overlap=$(wc -l < 06_peaks/qc_plots/rep1_overlapping_rep2.bed)
rep2_overlap=$(wc -l < 06_peaks/qc_plots/rep2_overlapping_rep1.bed)

echo "Replicate overlap analysis:"
echo "Rep1 total peaks: $rep1_total"
echo "Rep2 total peaks: $rep2_total"
echo "Rep1 peaks overlapping Rep2: $rep1_overlap ($(echo "scale=1; $rep1_overlap*100/$rep1_total" | bc)%)"
echo "Rep2 peaks overlapping Rep1: $rep2_overlap ($(echo "scale=1; $rep2_overlap*100/$rep2_total" | bc)%)"

# Good overlap: >50% for TF ChIP-seq
```

## Understanding Peak Calling Results

### Expected Results for ARF27:
- **Peak number**: 500-5000 peaks (depends on experimental conditions)
- **Peak length**: 100-500 bp (median ~200-300 bp)
- **High-confidence peaks**: 20-50% with score >100
- **Replicate overlap**: >50% for good quality data

### Peak Quality Indicators:

#### Good Quality Peaks:
- **Narrow, sharp peaks**: TF binding is specific
- **Consistent between replicates**: >50% overlap
- **Reasonable distribution**: Spread across all chromosomes
- **Range of scores**: Mix of high and moderate confidence peaks

#### Warning Signs:
- **Very broad peaks**: May indicate poor specificity
- **Low replicate overlap**: <30% suggests technical issues
- **Too many peaks**: >10,000 may indicate non-specific binding
- **Too few peaks**: <100 may indicate poor ChIP efficiency

### File Format Understanding:

#### narrowPeak format (BED6+4):
1. **Chromosome**: chr1, chr2, etc.
2. **Start**: Peak start coordinate (0-based)
3. **End**: Peak end coordinate (1-based)
4. **Name**: Peak name
5. **Score**: MACS2 score (-10*log10(qvalue))
6. **Strand**: Always "."
7. **signalValue**: Fold-change at peak summit
8. **pValue**: -log10(pvalue)
9. **qValue**: -log10(qvalue)
10. **peak**: Relative summit position

## Common Issues and Troubleshooting

### Issue 1: No Peaks Called

**Possible Causes:**
- Poor ChIP efficiency
- Inappropriate parameters
- Insufficient sequencing depth
- Wrong control file

**Diagnosis:**
```bash
# Check if any enrichment exists
macs2 callpeak \
    --treatment ARF27_ChIP_merged.bam \
    --control Input_control_merged.bam \
    --qvalue 0.1 \
    --nomodel --extsize 150 \
    --name ARF27_test_relaxed

# Check log files for errors
tail -50 logs/ARF27_ChIP_merged_peak_calling.log
```

**Solutions:**
- Relax FDR threshold (--qvalue 0.1)
- Check input files with `samtools view`
- Verify control file is appropriate

### Issue 2: Too Many Peaks (>10,000)

**Possible Causes:**
- Non-specific binding
- Contamination
- Too lenient parameters

**Solutions:**
```bash
# Use stricter parameters
macs2 callpeak \
    --qvalue 0.01 \
    --pvalue 1e-6 \
    # ... other parameters

# Filter peaks by score
awk '$5 > 100' peaks.narrowPeak > filtered_peaks.narrowPeak
```

### Issue 3: Very Broad Peaks

**Possible Causes:**
- Wrong parameters for TF ChIP-seq
- Histone contamination
- Poor antibody specificity

**Solutions:**
```bash
# Use TF-specific parameters
macs2 callpeak \
    --nomodel \
    --extsize 147 \
    --shift -73 \
    --call-summits
```

### Issue 4: Poor Replicate Overlap

**Diagnosis:**
```bash
# Check individual replicate quality
wc -l 06_peaks/individual_replicates/*_peaks.narrowPeak

# Look for batch effects
# Compare peak score distributions between replicates
```

**Solutions:**
- Check if one replicate has technical issues
- Use more lenient parameters for overlap analysis
- Consider IDR (Irreproducible Discovery Rate) analysis

### Issue 5: MACS2 Memory Issues

**Error**: Memory allocation error

**Solutions:**
```bash
# Reduce memory usage
macs2 callpeak --buffer-size 10000 \
    # ... other parameters

# Or split by chromosome
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
do
    macs2 callpeak \
        --treatment <(samtools view -b input.bam $chr) \
        --control <(samtools view -b control.bam $chr) \
        # ... other parameters
done
```

## Peak Calling Parameter Optimization

### For Different Analysis Goals:

#### High Sensitivity (find more peaks):
```bash
macs2 callpeak \
    --qvalue 0.1 \
    --nomodel \
    --extsize 150 \
    --shift -75
```

#### High Specificity (fewer, more confident peaks):
```bash
macs2 callpeak \
    --qvalue 0.01 \
    --pvalue 1e-6 \
    --nomodel \
    --extsize 150 \
    --shift -75
```

#### Broad region calling (if needed):
```bash
macs2 callpeak \
    --broad \
    --broad-cutoff 0.1 \
    --nomodel
```

## Analysis Log Update

```bash
# Document peak calling results
cat >> logs/step6_peak_calling_log.txt << EOF

Step 6: Peak Calling Results  
Date: $(date)

MACS2 Parameters Used:
- FDR threshold (qvalue): 0.05
- Genome size: 2.1e9 bp
- Extension size: 150 bp
- Model: No shifting model (--nomodel)
- Summits: Called precise summits

Peak Calling Summary:
- ARF27_ChIP_rep1: $(wc -l < 06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak) peaks
- ARF27_ChIP_rep2: $(wc -l < 06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak) peaks  
- ARF27_ChIP_merged: $(wc -l < 06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak) peaks

Quality Metrics:
- Replicate overlap: $(echo "scale=1; $(wc -l < 06_peaks/qc_plots/rep1_overlapping_rep2.bed)*100/$(wc -l < 06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak)" | bc)%
- High-confidence peaks (score>100): $(awk '$5 > 100' 06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak | wc -l)
- Median peak length: [From R analysis] bp

Files Generated:
- Individual replicate peaks: 2 files
- Merged replicate peaks: 1 file  
- Peak summits: 3 files
- QC plots and statistics: Multiple files

Issues Encountered: [None/List problems and solutions]

Quality Assessment: [Good/Acceptable/Poor] - [Brief justification]

Ready for Step 7: Peak Blacklist Filtering

EOF
```

## Key Takeaways

- MACS2 is the gold standard for ChIP-seq peak calling
- Transcription factor parameters differ from histone mark analysis
- Peak calling quality depends on experimental design and execution
- Replicate consistency is crucial for reliable results
- Parameter optimization may be needed for specific experimental conditions

## Next Steps Checklist

Before proceeding to Step 7 (Peak Blacklist Filtering):

- [ ] Peaks successfully called on individual replicates
- [ ] Peaks called on merged replicates  
- [ ] Peak statistics calculated and reasonable
- [ ] Replicate overlap assessed (>50% expected)
- [ ] Peak length distribution analyzed
- [ ] Output files verified and organized
- [ ] Any quality issues documented

Ready for **Step 7: Peak Blacklist Filtering**? We'll remove peaks that fall in problematic genomic regions to improve the specificity and reliability of our ARF27 binding site predictions.
