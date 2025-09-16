# Step 5: Data Normalization for ChIP-seq Analysis

## Objective
Normalize ChIP-seq data to enable quantitative comparisons between samples, replicates, and experimental conditions by accounting for differences in library size, sequencing depth, and technical variation.

## Why Normalization is Critical for ChIP-seq

### The Normalization Problem:
- **Different sequencing depths**: Samples may have 20M vs 50M reads
- **Library preparation variation**: Different PCR amplification efficiency
- **Sequencing bias**: Lane effects and batch differences
- **Biological variation**: Different amounts of chromatin input

### Impact on ARF27 Analysis:
- **Peak calling**: Unnormalized data leads to biased peak detection
- **Differential binding**: Cannot compare binding between conditions
- **Quantitative analysis**: Peak intensities not directly comparable
- **Visualization**: Browser tracks need consistent scaling

### Types of Normalization:

#### 1. Standard Library Size Normalization:
- **RPKM/FPKM**: Reads Per Kilobase per Million mapped reads
- **CPM**: Counts Per Million mapped reads
- **RPM**: Reads Per Million mapped reads
- **Library size scaling**: Scale to common read depth

#### 2. Spike-in Normalization:
- **External spike-ins**: Known amounts of non-endogenous DNA
- **Cross-species spike-ins**: E.g., Drosophila chromatin in plant samples
- **Accounts for**: Global changes in protein binding

## Understanding Normalization Methods

### Standard Normalization (Most Common):
```
Normalized_reads = (Raw_reads × 1,000,000) / Total_mapped_reads
```

### Spike-in Normalization (Advanced):
```
Scaling_factor = Spike_in_reads_sample / Mean_spike_in_reads_all_samples
Normalized_reads = Raw_reads / Scaling_factor
```

### When to Use Each Method:

#### Standard Normalization:
- **Most ChIP-seq experiments**: Comparing between samples/conditions
- **Assumption**: Total binding capacity is similar across samples
- **Good for**: ARF27 binding changes in specific conditions

#### Spike-in Normalization:
- **Global protein level changes**: When ARF27 total protein varies
- **Cross-species comparisons**: Different organisms
- **Temporal studies**: Development or stress responses
- **Requires**: Spike-in controls added during ChIP

## Step-by-Step Instructions

### 1. Install Required Tools

```bash
# Install deepTools for normalization and visualization
conda install -c bioconda deeptools

# Install R packages for statistical normalization
# (We'll use this later for differential analysis)
R -e "install.packages('BiocManager')"
R -e "BiocManager::install(c('DiffBind', 'csaw', 'edgeR'))"

# Verify installations
multiBamSummary --version
bamCoverage --version
```

### 2. Assess Library Sizes

First, examine the sequencing depth of each sample:

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Count reads in each deduplicated BAM file
echo "Sample\tTotal_Reads\tMillion_Reads" > logs/library_sizes.txt

for bam in 05_processed_bams/deduplicated/*_dedup.bam
do
    sample=$(basename $bam _dedup.bam)
    total_reads=$(samtools view -c $bam)
    million_reads=$(echo "scale=2; $total_reads / 1000000" | bc)
    echo -e "${sample}\t${total_reads}\t${million_reads}" >> logs/library_sizes.txt
done

# Display library sizes
echo "Library sizes (millions of reads):"
column -t logs/library_sizes.txt
```

### 3. Create Normalized BigWig Files

BigWig files are essential for visualization and quantitative analysis:

#### Standard RPM Normalization:

```bash
# Create directory for normalized data
mkdir -p 05_processed_bams/normalized_bigwigs

# Generate RPM-normalized BigWig files for each sample
echo "Creating RPM-normalized BigWig files..."

# ARF27 ChIP replicate 1
bamCoverage \
    --bam 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_rep1_RPM.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150

# ARF27 ChIP replicate 2
bamCoverage \
    --bam 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_rep2_RPM.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150

# Input control replicate 1
bamCoverage \
    --bam 05_processed_bams/deduplicated/Input_control_rep1_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/Input_control_rep1_RPM.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150

# Input control replicate 2
bamCoverage \
    --bam 05_processed_bams/deduplicated/Input_control_rep2_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/Input_control_rep2_RPM.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150
```

**Command Explanation:**
- `--bam`: Input deduplicated BAM file
- `--outFileName`: Output BigWig file name
- `--binSize 10`: 10bp resolution (good for TF ChIP-seq)
- `--normalizeUsing RPM`: Reads Per Million normalization
- `--effectiveGenomeSize`: Maize genome size (~2.1 Gb)
- `--numberOfProcessors 8`: Use 8 CPU cores
- `--extendReads 150`: Extend reads to fragment length

### 4. Alternative: Create Input-Subtracted Tracks

For better visualization of ChIP enrichment:

```bash
# Create input-subtracted normalized tracks
echo "Creating input-subtracted BigWig files..."

# ARF27 ChIP rep1 minus Input rep1
bamCompare \
    --bamfile1 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    --bamfile2 05_processed_bams/deduplicated/Input_control_rep1_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_rep1_vs_Input_log2.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --operation log2 \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150

# ARF27 ChIP rep2 minus Input rep2
bamCompare \
    --bamfile1 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam \
    --bamfile2 05_processed_bams/deduplicated/Input_control_rep2_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_rep2_vs_Input_log2.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --operation log2 \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8 \
    --extendReads 150
```

### 5. Spike-in Normalization (If Applicable)

If you have spike-in controls (e.g., Drosophila chromatin):

```bash
# Map reads to spike-in genome (example with Drosophila)
# First, create combined reference (maize + spike-in)
# cat maize_B73_v4.fa drosophila_dm6.fa > combined_reference.fa
# bwa index combined_reference.fa

# Count spike-in reads in each sample
echo "Sample\tSpike_in_Reads\tEndogenous_Reads\tSpike_in_Percent" > logs/spike_in_counts.txt

for bam in 05_processed_bams/deduplicated/*_dedup.bam
do
    sample=$(basename $bam _dedup.bam)
    
    # Count reads mapping to spike-in chromosomes (dm6 chromosomes start with "chr")
    spike_reads=$(samtools view -c $bam chr2L chr2R chr3L chr3R chr4 chrX chrY 2>/dev/null || echo "0")
    
    # Count total reads
    total_reads=$(samtools view -c $bam)
    endogenous_reads=$((total_reads - spike_reads))
    
    # Calculate percentage
    spike_percent=$(echo "scale=2; $spike_reads * 100 / $total_reads" | bc)
    
    echo -e "${sample}\t${spike_reads}\t${endogenous_reads}\t${spike_percent}" >> logs/spike_in_counts.txt
done

# Calculate normalization factors
R --vanilla << 'EOF'
# Read spike-in data
spike_data <- read.table("logs/spike_in_counts.txt", header=TRUE, sep="\t")

# Calculate normalization factors
mean_spike <- mean(spike_data$Spike_in_Reads)
spike_data$Norm_Factor <- spike_data$Spike_in_Reads / mean_spike

# Save normalization factors
write.table(spike_data[,c("Sample", "Norm_Factor")], 
           "logs/spike_in_norm_factors.txt", 
           sep="\t", quote=FALSE, row.names=FALSE)

print("Spike-in normalization factors:")
print(spike_data[,c("Sample", "Spike_in_Reads", "Norm_Factor")])
EOF

# Apply spike-in normalization to create BigWig files
# (This would require custom scaling - shown conceptually)
```

### 6. Quality Control of Normalized Data

#### A. Compare Library Sizes After Normalization:

```bash
# Create comparison plot of library sizes
R --vanilla << 'EOF'
# Read library size data
lib_sizes <- read.table("logs/library_sizes.txt", header=TRUE, sep="\t")

# Create simple bar plot
png("logs/library_sizes_comparison.png", width=800, height=600)
barplot(lib_sizes$Million_Reads, 
        names.arg=lib_sizes$Sample, 
        las=2, 
        main="Library Sizes (Million Reads)",
        ylab="Million Reads",
        col=c("lightblue", "lightblue", "lightcoral", "lightcoral"))
dev.off()

# Print summary
print("Library size summary:")
print(summary(lib_sizes$Million_Reads))
EOF
```

#### B. Check BigWig File Integrity:

```bash
# Verify BigWig files were created successfully
echo "Checking BigWig file integrity..."
for bw in 05_processed_bams/normalized_bigwigs/*.bw
do
    if [ -f "$bw" ]; then
        size=$(ls -lh "$bw" | awk '{print $5}')
        echo "✓ $(basename $bw): $size"
    else
        echo "✗ $(basename $bw): Missing"
    fi
done
```

#### C. Sample Correlation Analysis:

```bash
# Create correlation matrix of normalized samples
echo "Creating sample correlation matrix..."

multiBamSummary bins \
    --bamfiles 05_processed_bams/deduplicated/*_dedup.bam \
    --labels ARF27_rep1 ARF27_rep2 Input_rep1 Input_rep2 \
    --numberOfProcessors 8 \
    --binSize 10000 \
    --outFileName logs/sample_correlation_matrix.npz

# Create correlation heatmap
plotCorrelation \
    --corData logs/sample_correlation_matrix.npz \
    --plotFile logs/sample_correlation_heatmap.png \
    --corMethod pearson \
    --whatToPlot heatmap \
    --plotNumbers \
    --plotTitle "Sample Correlation Matrix"

# Create PCA plot
plotPCA \
    --corData logs/sample_correlation_matrix.npz \
    --plotFile logs/sample_PCA.png \
    --plotTitle "Principal Component Analysis"
```

## Interpreting Normalization Results

### Expected Correlation Patterns:

#### Good Quality Data:
- **Biological replicates**: r > 0.8 (ChIP samples)
- **Input replicates**: r > 0.9 (less specific)
- **ChIP vs Input**: r = 0.3-0.7 (some correlation expected)

#### Warning Signs:
- **Low replicate correlation**: r < 0.7 (technical issues)
- **High ChIP-Input correlation**: r > 0.8 (poor enrichment)
- **Batch effects**: Samples cluster by prep date, not biology

### Library Size Assessment:

```bash
# Calculate coefficient of variation for library sizes
R --vanilla << 'EOF'
lib_sizes <- read.table("logs/library_sizes.txt", header=TRUE, sep="\t")

# Calculate CV for ChIP samples
chip_sizes <- lib_sizes[grep("ARF27_ChIP", lib_sizes$Sample), "Million_Reads"]
chip_cv <- sd(chip_sizes) / mean(chip_sizes) * 100

# Calculate CV for Input samples  
input_sizes <- lib_sizes[grep("Input", lib_sizes$Sample), "Million_Reads"]
input_cv <- sd(input_sizes) / mean(input_sizes) * 100

cat("ChIP library size CV:", round(chip_cv, 1), "%\n")
cat("Input library size CV:", round(input_cv, 1), "%\n")

# Good: CV < 20%
# Acceptable: CV < 50%
# Poor: CV > 50%
EOF
```

## Advanced Normalization Options

### 1. GC Content Correction:

If you notice GC bias in your data:

```bash
# Create GC-corrected BigWig files
bamCoverage \
    --bam 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam \
    --outFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_rep1_GC_corrected.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --normalizeUsing RPM \
    --correctGCBias \
    --genome ../reference_genome/maize_B73_v4.fa \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8
```

### 2. Blacklist Region Exclusion:

Exclude problematic genomic regions from normalization:

```bash
# If you have a blacklist BED file (we'll create one in Step 7)
# bamCoverage \
#     --bam input.bam \
#     --outFileName output.bw \
#     --blackListFileName maize_blacklist.bed \
#     --normalizeUsing RPM
```

### 3. Bin Size Optimization:

For different analysis purposes:

```bash
# High resolution for detailed analysis (slower)
bamCoverage --binSize 1 ...   # 1bp resolution

# Medium resolution for standard analysis
bamCoverage --binSize 10 ...  # 10bp resolution (recommended)

# Low resolution for genome-wide overview
bamCoverage --binSize 100 ... # 100bp resolution (faster)
```

## Common Issues and Troubleshooting

### Issue 1: Large Differences in Library Sizes

**Problem**: Some samples have 2x+ more reads than others

**Assessment**:
```bash
# Check the range of library sizes
awk 'NR>1 {print $3}' logs/library_sizes.txt | sort -n
```

**Solutions**:
- **If systematic**: May indicate batch effects - document in methods
- **If random**: Normal biological/technical variation - normalization handles this
- **If extreme** (>5x difference): Consider subsampling largest libraries

### Issue 2: Poor Replicate Correlation

**Problem**: Biological replicates correlate poorly (r < 0.7)

**Diagnosis**:
```bash
# Check correlation values
grep -A 10 "Correlation matrix" logs/sample_correlation_matrix.txt
```

**Solutions**:
- Check for batch effects in sample preparation
- Examine individual sample quality metrics
- Consider excluding poor-quality replicates
- Investigate experimental conditions

### Issue 3: BigWig Generation Fails

**Common errors**:
- Memory issues with large BAM files
- Incorrect effective genome size
- Missing chromosome names

**Solutions**:
```bash
# Reduce memory usage
bamCoverage --binSize 50 ...  # Larger bins use less memory

# Check chromosome names match
samtools view -H 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam | grep "@SQ"

# For memory issues
export TMPDIR=/path/to/large/disk/tmp
```

### Issue 4: Spike-in Controls Show High Variation

**Problem**: Spike-in percentages vary >2-fold between samples

**Possible causes**:
- Inconsistent spike-in addition
- Cross-contamination
- Processing errors

**Solutions**:
- Document variation and proceed with caution
- Consider excluding samples with extreme spike-in values
- Use standard normalization as backup

## File Organization and Quality Checks

### 1. Verify Output Files:

```bash
# Check normalized data directory
ls -la 05_processed_bams/normalized_bigwigs/

# Expected files:
# - *_RPM.bw (RPM normalized tracks)
# - *_vs_Input_log2.bw (input-subtracted tracks)
# - Quality control plots in logs/
```

### 2. File Size Expectations:

```bash
# BigWig files should be much smaller than BAM files
echo "File size comparison:"
echo "BAM files:"
ls -lh 05_processed_bams/deduplicated/*.bam | awk '{print $5, $9}'
echo "BigWig files:"
ls -lh 05_processed_bams/normalized_bigwigs/*.bw | awk '{print $5, $9}'
```

**Expected sizes**:
- **BAM files**: 5-20 GB
- **BigWig files**: 50-500 MB (10-50x smaller)

## Analysis Log Update

```bash
# Document normalization results
cat >> logs/step5_normalization_log.txt << EOF

Step 5: Data Normalization Results
Date: $(date)

Library Sizes (Million Reads):
$(cat logs/library_sizes.txt)

Normalization Method: RPM (Reads Per Million)
Effective Genome Size: 2.1 Gb (maize B73-v4)
Bin Size: 10 bp
Read Extension: 150 bp

Quality Control:
- BigWig files generated successfully: [Yes/No]
- Sample correlations calculated: [Yes/No]
- Library size variation acceptable: [Yes/No]
- Spike-in normalization applied: [Yes/No/N/A]

Files Generated:
- RPM-normalized BigWig tracks: $(ls 05_processed_bams/normalized_bigwigs/*_RPM.bw | wc -l) files
- Input-subtracted tracks: $(ls 05_processed_bams/normalized_bigwigs/*_vs_Input_log2.bw | wc -l) files
- Correlation analysis: Complete

Issues Encountered: [None/List problems and solutions]

Ready for Step 6: Peak Calling

EOF
```

## Key Takeaways

- Normalization enables quantitative comparison between samples
- RPM normalization is standard for most ChIP-seq analyses
- BigWig files are essential for visualization and downstream analysis
- Sample correlation analysis reveals data quality and experimental issues
- Proper normalization is critical for accurate peak calling and differential analysis

## Next Steps Checklist

Before proceeding to Step 6 (Peak Calling):

- [ ] All samples successfully normalized to RPM
- [ ] BigWig files generated and verified
- [ ] Sample correlation analysis completed
- [ ] Library size variation documented
- [ ] Any quality issues identified and addressed
- [ ] Input-subtracted tracks created for visualization

Ready for **Step 6: Peak Calling**? We'll use MACS2 to identify ARF27 binding sites by calling peaks from the normalized ChIP-seq data, using appropriate parameters optimized for transcription factor ChIP-seq.
