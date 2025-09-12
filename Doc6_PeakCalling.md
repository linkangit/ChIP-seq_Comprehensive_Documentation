# Document 6: Peak Calling with MACS2

## Understanding Peak Calling in ChIP-seq

Peak calling is the process of identifying genomic regions where your protein of interest (ARF27) is significantly enriched compared to background levels. Think of it as finding the "hotspots" where ARF27 binds to DNA.

After alignment, you have millions of reads scattered across the genome. Most of these reads represent background binding or random sampling, but some genomic regions will have many more reads than expected by chance. These enriched regions are called "peaks" and represent likely binding sites for ARF27.

### The Biology Behind ChIP-seq Peaks

In your ARF27 experiment, a true binding site creates a characteristic pattern:

**Signal Enrichment**: Where ARF27 binds, you'll see more sequencing reads than in random genomic regions. This enrichment occurs because:
- Cross-linking captures ARF27-DNA complexes
- Immunoprecipitation specifically pulls down ARF27-bound DNA fragments
- These fragments are preferentially sequenced

**Peak Shape**: True transcription factor binding sites typically show:
- A distinct summit (highest point) indicating the most likely binding location
- Gradual decrease in read density moving away from the binding site
- Peak width typically 100-500 base pairs for transcription factors

**Background vs. Signal**: Not all read accumulations represent true binding:
- Random sampling creates some read clusters by chance
- Repetitive sequences may accumulate reads artificially
- Technical artifacts can create false peaks

### Why We Need Statistical Analysis

Simply counting reads isn't sufficient because:
- Random variation creates apparent "peaks" by chance
- Different genomic regions have different baseline read densities
- Sequencing depth varies between samples
- We need to distinguish true biological signal from technical noise

MACS2 (Model-based Analysis of ChIP-Seq) solves these problems using sophisticated statistical models.

## How MACS2 Works

MACS2 uses several steps to identify genuine binding sites:

### 1. Background Modeling
MACS2 estimates the expected background read density across the genome by:
- Analyzing read distribution in control samples (input DNA)
- Modeling local sequence bias and mappability
- Accounting for global sequencing depth differences

### 2. Fragment Size Estimation
For single-end sequencing, MACS2 estimates the original DNA fragment size by:
- Finding pairs of forward and reverse reads
- Calculating the distance between strand-specific read clusters
- Using this to extend reads to full fragment length

### 3. Peak Detection
MACS2 identifies candidate peaks by:
- Scanning the genome for regions with elevated read density
- Comparing ChIP signal to control background
- Calculating statistical significance for each potential peak

### 4. Statistical Testing
For each candidate peak, MACS2:
- Calculates a p-value based on background expectation
- Applies multiple testing correction (false discovery rate)
- Reports q-values (FDR-corrected p-values)

## Setting Up Peak Calling

```bash
# Make sure we're in the right directory and environment
cd chipseq_analysis
conda activate chipseq

# Create directory for peak calling results
mkdir -p peaks
```

## Step 1: Identify Your Sample Types

Before peak calling, organize your samples by type:

```bash
# Look at your aligned BAM files
ls alignment/*_filtered.bam

# Identify ChIP and control samples
echo "ChIP samples (should contain your protein of interest):"
ls alignment/*ChIP*_filtered.bam 2>/dev/null || echo "No files matching *ChIP* pattern"
ls alignment/*ARF27*_filtered.bam 2>/dev/null || echo "No files matching *ARF27* pattern"

echo "Control samples (input DNA or IgG controls):"
ls alignment/*input*_filtered.bam 2>/dev/null || echo "No files matching *input* pattern"
ls alignment/*IgG*_filtered.bam 2>/dev/null || echo "No files matching *IgG* pattern"
```

## Step 2: Basic Peak Calling (Single Sample)

Let's start with a simple example using one ChIP sample against one control:

```bash
# Example: Call peaks for sample1 ChIP vs input
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.05 \
    --keep-dup all \
    --call-summits
```

### Understanding MACS2 Parameters

**-t alignment/sample1_ChIP_filtered.bam**
- Treatment file (your ChIP sample)
- Contains reads from ARF27-bound DNA fragments
- Should show enrichment at true binding sites

**-c alignment/sample1_input_filtered.bam**
- Control file (background DNA)
- Represents random genomic sampling
- Used to model background read distribution
- Essential for distinguishing real signal from noise

**-n sample1_ChIP**
- Name prefix for output files
- MACS2 will create files like: sample1_ChIP_peaks.narrowPeak
- Use descriptive names to track your samples

**--outdir peaks**
- Directory where output files will be created
- Keeps results organized and separate from input data

**-g 2.1e9**
- Mappable genome size for maize (2.1 billion base pairs)
- Used for background rate calculations
- Critical for accurate statistical testing
- Slightly smaller than total genome size (2.3e9) due to repetitive regions

**--nomodel**
- Skip fragment size modeling step
- Appropriate for transcription factor ChIP-seq
- Use when you have reliable fragment size estimates
- Speeds up analysis and reduces artifacts

**--extsize 147**
- Extend reads to this fragment size (147 base pairs)
- Based on nucleosome protection length
- Reasonable default for transcription factor ChIP-seq
- Can adjust based on your experimental fragment size

**-q 0.05**
- Q-value (FDR) threshold for peak significance
- 0.05 = 5% false discovery rate
- More stringent than p-value due to multiple testing correction
- Standard threshold balancing sensitivity and specificity

**--keep-dup all**
- Keep all duplicate reads
- Important for ChIP-seq where true enrichment creates "duplicates"
- Alternative: --keep-dup auto (let MACS2 decide)

**--call-summits**
- Identify the most likely binding position within each peak
- Creates _summits.bed file with point locations
- Useful for motif analysis and precise binding site identification

## Step 3: Call Peaks for Multiple Samples

Process each ChIP sample individually:

```bash
# Call peaks for sample 2
macs2 callpeak \
    -t alignment/sample2_ChIP_filtered.bam \
    -c alignment/sample2_input_filtered.bam \
    -n sample2_ChIP \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.05 \
    --keep-dup all \
    --call-summits

# Call peaks for sample 3  
macs2 callpeak \
    -t alignment/sample3_ChIP_filtered.bam \
    -c alignment/sample3_input_filtered.bam \
    -n sample3_ChIP \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.05 \
    --keep-dup all \
    --call-summits

# Continue for all your ChIP samples...
```

## Step 4: Understanding MACS2 Output Files

MACS2 creates several output files for each sample:

```bash
# Look at the output files for sample1
ls peaks/sample1_ChIP*

# Examine the different file types
echo "File types created by MACS2:"
ls peaks/*.narrowPeak 2>/dev/null && echo "- narrowPeak: Main peak coordinates"
ls peaks/*.summits.bed 2>/dev/null && echo "- summits.bed: Peak summit positions"  
ls peaks/*.xls 2>/dev/null && echo "- xls: Detailed peak information"
ls peaks/*.bdg 2>/dev/null && echo "- bdg: Bedgraph coverage tracks"
```

### Key Output Files Explained

**sample1_ChIP_peaks.narrowPeak**
- Main output file with peak coordinates
- Standard BED format with additional ChIP-seq specific columns
- This is the file you'll use for most downstream analyses

**sample1_ChIP_summits.bed**
- Point locations of peak summits
- Represents the most likely binding positions
- Useful for motif discovery and precise annotation

**sample1_ChIP_peaks.xls** 
- Human-readable spreadsheet with detailed peak information
- Includes p-values, q-values, fold enrichment
- Good for manual inspection and filtering

### Examining Peak Files

```bash
# Look at the first few peaks
head -5 peaks/sample1_ChIP_peaks.narrowPeak

# Count how many peaks were called
wc -l peaks/sample1_ChIP_peaks.narrowPeak

# Look at peak summits
head -5 peaks/sample1_ChIP_summits.bed
```

### Understanding the narrowPeak Format

The narrowPeak format contains 10 columns:
1. **Chromosome**: Which chromosome the peak is on
2. **Start**: Peak start coordinate  
3. **End**: Peak end coordinate
4. **Name**: Peak identifier
5. **Score**: Integer score (0-1000)
6. **Strand**: Strand information (usually ".")
7. **Signal Value**: Measurement of peak strength
8. **P-value**: Statistical significance (-log10 p-value)
9. **Q-value**: FDR-corrected significance (-log10 q-value)  
10. **Peak**: Relative summit position within the peak

## Step 5: Evaluating Peak Calling Results

### Basic Peak Statistics

```bash
# Count peaks for each sample
echo "Peak counts per sample:"
for peak_file in peaks/*.narrowPeak; do
    sample=$(basename $peak_file _peaks.narrowPeak)
    count=$(wc -l < $peak_file)
    echo "$sample: $count peaks"
done
```

### Peak Quality Assessment

```bash
# Look at peak score distribution for sample1
cut -f7 peaks/sample1_ChIP_peaks.narrowPeak | sort -n | tail -20

# Look at q-value distribution (higher values = more significant)
cut -f9 peaks/sample1_ChIP_peaks.narrowPeak | sort -n | tail -20

# Count highly significant peaks (q-value > 20, meaning q < 1e-20)
awk '$9 > 20' peaks/sample1_ChIP_peaks.narrowPeak | wc -l
```

### Expected Results for Good ChIP-seq

**Peak Numbers**:
- **1,000-5,000 peaks**: Conservative, high-confidence binding sites
- **5,000-20,000 peaks**: Moderate stringency, good coverage
- **20,000+ peaks**: Permissive, may include weaker sites

**Peak Quality Indicators**:
- Wide range of peak scores (some very high values)
- Many peaks with q-values > 10 (q < 1e-10)
- Reasonable peak widths (100-1000 bp for transcription factors)

## Step 6: More Stringent Peak Calling

If you get too many peaks, try more stringent parameters:

```bash
# More stringent peak calling (lower FDR threshold)
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_stringent \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.01 \
    --keep-dup all \
    --call-summits

# Count peaks with stringent parameters
wc -l peaks/sample1_ChIP_stringent_peaks.narrowPeak
```

## Step 7: Peak Calling Without Control (if necessary)

If you don't have control samples, MACS2 can still call peaks:

```bash
# Peak calling without control (not ideal but sometimes necessary)
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -n sample1_ChIP_nocontrol \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.05 \
    --keep-dup all \
    --call-summits

# Note: This will likely give more false positives
```

## Step 8: Merging Biological Replicates

If you have multiple biological replicates, you can merge them for increased power:

```bash
# Merge BAM files from biological replicates
samtools merge alignment/merged_ChIP.bam \
    alignment/sample1_ChIP_filtered.bam \
    alignment/sample2_ChIP_filtered.bam

samtools merge alignment/merged_input.bam \
    alignment/sample1_input_filtered.bam \
    alignment/sample2_input_filtered.bam

# Index merged files
samtools index alignment/merged_ChIP.bam
samtools index alignment/merged_input.bam

# Call peaks on merged data
macs2 callpeak \
    -t alignment/merged_ChIP.bam \
    -c alignment/merged_input.bam \
    -n merged_replicates \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.05 \
    --keep-dup all \
    --call-summits
```

## Step 9: Alternative Peak Calling Parameters

### For Broader Peaks (if needed)

```bash
# Call broader peaks (for some transcription factors)
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_broad \
    --outdir peaks \
    -g 2.1e9 \
    --broad \
    --broad-cutoff 0.1 \
    --keep-dup all
```

### Fragment Size Optimization

```bash
# Let MACS2 estimate fragment size automatically
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_auto \
    --outdir peaks \
    -g 2.1e9 \
    -q 0.05 \
    --keep-dup all \
    --call-summits

# Use custom fragment size if you know it
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_custom \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 200 \
    -q 0.05 \
    --keep-dup all \
    --call-summits
```

## Step 10: Quality Control and Validation

### Create Peak Summary

```bash
# Create comprehensive peak summary
cat > peaks/peak_calling_summary.txt << 'EOF'
PEAK CALLING SUMMARY
===================

MACS2 Parameters Used:
- Genome size: 2.1e9 (maize mappable genome)
- Q-value threshold: 0.05 (5% FDR)
- Fragment size: 147 bp (--extsize)
- Model: --nomodel (skip fragment modeling)
- Duplicates: kept all

Peak Counts by Sample:
EOF

# Add peak counts to summary
for peak_file in peaks/*.narrowPeak; do
    sample=$(basename $peak_file _peaks.narrowPeak)
    count=$(wc -l < $peak_file)
    echo "$sample: $count peaks" >> peaks/peak_calling_summary.txt
done

# View the summary
cat peaks/peak_calling_summary.txt
```

### Validate Peak Quality

```bash
# Check for very large peaks (may indicate artifacts)
echo "Checking for unusually large peaks:"
for peak_file in peaks/*.narrowPeak; do
    sample=$(basename $peak_file _peaks.narrowPeak)
    large_peaks=$(awk '$3-$2 > 5000' $peak_file | wc -l)
    echo "$sample: $large_peaks peaks larger than 5kb"
done

# Check peak width distribution
echo "Peak width statistics for sample1:"
awk '{print $3-$2}' peaks/sample1_ChIP_peaks.narrowPeak | sort -n | awk '
BEGIN { count = 0; sum = 0 }
{ 
    widths[count] = $1; 
    sum += $1; 
    count++ 
}
END { 
    print "Count:", count
    print "Mean width:", sum/count
    print "Median width:", widths[int(count/2)]
    print "Min width:", widths[0]
    print "Max width:", widths[count-1]
}'
```

## Troubleshooting Peak Calling Issues

### Too Few Peaks

**Possible causes and solutions**:

```bash
# Try less stringent q-value
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_relaxed \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.1 \
    --keep-dup all \
    --call-summits

# Check if you have enough reads
samtools flagstat alignment/sample1_ChIP_filtered.bam
```

### Too Many Peaks

```bash
# Use more stringent threshold
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -c alignment/sample1_input_filtered.bam \
    -n sample1_ChIP_strict \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.01 \
    --keep-dup all \
    --call-summits

# Filter peaks by additional criteria
awk '$7 > 10 && $9 > 5' peaks/sample1_ChIP_peaks.narrowPeak > peaks/sample1_ChIP_filtered_peaks.bed
```

### No Control Sample

```bash
# If you must proceed without control, use these parameters
macs2 callpeak \
    -t alignment/sample1_ChIP_filtered.bam \
    -n sample1_ChIP_no_control \
    --outdir peaks \
    -g 2.1e9 \
    --nomodel \
    --extsize 147 \
    -q 0.001 \
    --keep-dup all \
    --call-summits

# Be extra stringent without control
```

## Understanding Peak Calling Success

### Good Peak Calling Results

**Peak Numbers**: 1,000-20,000 peaks (depending on stringency)
**Peak Quality**: Wide range of significance values
**Peak Sizes**: Mostly 100-1,000 bp for transcription factors
**Reproducibility**: Similar peak numbers between biological replicates

### Warning Signs

**Very few peaks** (<100): May indicate poor enrichment or overly stringent parameters
**Too many peaks** (>50,000): May indicate weak enrichment or contamination
**Very large peaks** (>10,000 bp): May indicate broad chromatin domains rather than specific binding
**No high-confidence peaks**: May indicate experimental failure

## Preparing for Peak Annotation

You now have peak files that represent likely ARF27 binding sites. The next step will be to annotate these peaks to understand:
- Which genes are near the binding sites
- Whether peaks are in promoters, gene bodies, or intergenic regions
- What biological processes might be regulated by ARF27

### Files Ready for Annotation

```bash
# Verify you have the peak files needed
ls peaks/*.narrowPeak
ls peaks/*.summits.bed

# Check that files are not empty
for file in peaks/*.narrowPeak; do
    lines=$(wc -l < $file)
    echo "$(basename $file): $lines peaks"
done
```

Your peak calling is now complete, and you have identified genomic regions where ARF27 likely binds to regulate gene expression in maize.
