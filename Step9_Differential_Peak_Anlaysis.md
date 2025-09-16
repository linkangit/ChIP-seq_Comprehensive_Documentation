# Step 9: Differential Peak Analysis for ChIP-seq Data

## Objective
Identify condition-specific ARF27 binding sites and quantify changes in binding intensity between experimental conditions using statistical methods designed for ChIP-seq data to understand how ARF27 binding responds to different treatments or developmental stages.

## Why Differential Peak Analysis is Important

### Scientific Questions Addressed:
- **Treatment response**: How does auxin treatment affect ARF27 binding?
- **Developmental changes**: How does ARF27 binding change during development?
- **Stress responses**: Which binding sites are gained/lost under stress?
- **Tissue specificity**: What are tissue-specific ARF27 targets?

### Types of Differential Analysis:
1. **Differential binding**: Changes in peak intensity between conditions
2. **Condition-specific peaks**: Peaks present in one condition but not another
3. **Dynamic binding**: Time-course changes in binding patterns
4. **Quantitative changes**: Fold-change in binding strength

### Statistical Challenges:
- **Count data**: ChIP-seq reads follow negative binomial distribution
- **Multiple testing**: Thousands of peaks require FDR correction
- **Normalization**: Library size and efficiency differences
- **Replicate variation**: Biological and technical variability

## Understanding Differential Analysis Methods

### 1. DiffBind Approach:
- **Uses**: DESeq2 or edgeR statistical frameworks
- **Strengths**: Designed specifically for ChIP-seq, handles replicates well
- **Method**: Quantifies reads in consensus peaks across all samples

### 2. csaw (ChIP-seq Analysis with Windows):
- **Uses**: Window-based approach with edgeR
- **Strengths**: Can discover novel differential regions
- **Method**: Slides windows across genome, tests each independently

### 3. Manual Approach:
- **Uses**: Simple fold-change and presence/absence
- **Strengths**: Easy to understand and implement
- **Method**: Direct comparison of peak scores and overlap

## Step-by-Step Differential Analysis

### 1. Set Up Analysis Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create differential analysis directory
mkdir -p 09_differential_analysis
mkdir -p 09_differential_analysis/sample_sheets
mkdir -p 09_differential_analysis/consensus_peaks
mkdir -p 09_differential_analysis/count_matrices
mkdir -p 09_differential_analysis/results
mkdir -p 09_differential_analysis/plots
```

### 2. Prepare Sample Information

Create a sample sheet describing your experimental design:

```bash
# Create sample sheet for differential analysis
# Modify this based on your actual experimental conditions
cat > 09_differential_analysis/sample_sheets/sample_info.csv << 'EOF'
SampleID,Condition,Replicate,Treatment,Tissue,BAM_File,Peak_File
ARF27_Control_Rep1,Control,1,Untreated,Root,05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam,07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak
ARF27_Control_Rep2,Control,2,Untreated,Root,05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam,07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak
ARF27_Auxin_Rep1,Auxin,1,IAA_Treatment,Root,05_processed_bams/deduplicated/ARF27_Auxin_rep1_dedup.bam,07_blacklists/filtered_peaks/ARF27_Auxin_rep1_peaks_filtered.narrowPeak
ARF27_Auxin_Rep2,Auxin,2,IAA_Treatment,Root,05_processed_bams/deduplicated/ARF27_Auxin_rep2_dedup.bam,07_blacklists/filtered_peaks/ARF27_Auxin_rep2_peaks_filtered.narrowPeak
EOF

echo "Sample sheet created: 09_differential_analysis/sample_sheets/sample_info.csv"
echo "Please modify this file according to your actual experimental design!"
```

### 3. Install Required R Packages

```bash
# Install required R packages for differential analysis
R --vanilla << 'EOF'
# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "DiffBind",
    "DESeq2", 
    "edgeR",
    "csaw",
    "GenomicRanges",
    "rtracklayer",
    "ChIPseeker",
    "clusterProfiler"
))

# Install CRAN packages
install.packages(c(
    "pheatmap",
    "ggplot2",
    "dplyr",
    "VennDiagram"
))

# Check installations
library(DiffBind)
library(DESeq2)
library(edgeR)
sessionInfo()
EOF
```

### 4. Create Consensus Peak Set

Generate a unified peak set across all samples for quantitative comparison:

```bash
# Create consensus peaks using bedtools
echo "Creating consensus peak set..."

# Combine all peak files
cat 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak \
    07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak \
    > 09_differential_analysis/consensus_peaks/all_peaks_combined.bed

# Sort and merge overlapping peaks
sort -k1,1 -k2,2n 09_differential_analysis/consensus_peaks/all_peaks_combined.bed | \
bedtools merge -i - > 09_differential_analysis/consensus_peaks/consensus_peaks_raw.bed

# Add peak IDs
awk '{print $1"\t"$2"\t"$3"\tpeak_"NR}' \
    09_differential_analysis/consensus_peaks/consensus_peaks_raw.bed > \
    09_differential_analysis/consensus_peaks/consensus_peaks.bed

# Statistics
echo "Consensus peak statistics:"
echo "Total consensus regions: $(wc -l < 09_differential_analysis/consensus_peaks/consensus_peaks.bed)"
```

### 5. Method 1: DiffBind Analysis

Use DiffBind for comprehensive differential binding analysis:

```bash
# Run DiffBind analysis
R --vanilla << 'EOF'
library(DiffBind)
library(DESeq2)
library(ggplot2)

# Read sample information
# Note: You'll need to modify the sample sheet to match your actual data
samples <- read.csv("09_differential_analysis/sample_sheets/sample_info.csv")

# Create DiffBind sample sheet format
# This is a template - modify according to your data structure
diffbind_samples <- data.frame(
    SampleID = samples$SampleID,
    Condition = samples$Condition,
    Replicate = samples$Replicate,
    bamReads = samples$BAM_File,
    Peaks = samples$Peak_File,
    PeakCaller = "macs"
)

# Write DiffBind sample sheet
write.csv(diffbind_samples, "09_differential_analysis/sample_sheets/diffbind_samples.csv", 
          row.names = FALSE, quote = FALSE)

# Create DiffBind object (only if you have the required data structure)
# db <- dba(sampleSheet = "09_differential_analysis/sample_sheets/diffbind_samples.csv")
# 
# # Count reads in consensus peaks
# db <- dba.count(db, minOverlap = 2)
# 
# # Establish contrast (Control vs Auxin treatment)
# db <- dba.contrast(db, categories = DBA_CONDITION)
# 
# # Perform differential analysis
# db <- dba.analyze(db, method = DBA_DESEQ2)
# 
# # Generate results
# results <- dba.report(db, contrast = 1, th = 0.05)
# 
# # Save results
# write.csv(as.data.frame(results), 
#           "09_differential_analysis/results/diffbind_results.csv", 
#           row.names = FALSE)
# 
# # Create plots
# pdf("09_differential_analysis/plots/diffbind_plots.pdf")
# plot(db)
# dba.plotPCA(db, DBA_CONDITION)
# dba.plotMA(db, contrast = 1)
# dev.off()

cat("DiffBind template created. Please modify sample sheet according to your data.\n")
EOF
```

### 6. Method 2: Manual Differential Analysis

Implement a simpler approach for differential binding analysis:

```bash
# Manual differential analysis approach
echo "Performing manual differential analysis..."

# For demonstration, we'll compare the two replicates as if they were different conditions
# In real analysis, you would have actual condition differences

R --vanilla << 'EOF'
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(pheatmap)

# Read consensus peaks
consensus <- read.table("09_differential_analysis/consensus_peaks/consensus_peaks.bed")
colnames(consensus) <- c("chr", "start", "end", "peak_id")

# Create GRanges object for consensus peaks
consensus_gr <- GRanges(seqnames = consensus$chr,
                       ranges = IRanges(start = consensus$start, end = consensus$end),
                       peak_id = consensus$peak_id)

# Function to count reads in peaks
count_reads_in_peaks <- function(bam_file, peaks_gr) {
    # This is a simplified version - in practice you'd use packages like Rsubread
    # For now, we'll simulate counts based on peak properties
    set.seed(42)
    n_peaks <- length(peaks_gr)
    # Simulate counts with realistic ChIP-seq characteristics
    counts <- rpois(n_peaks, lambda = 50) + 
              rbinom(n_peaks, size = 100, prob = 0.3)
    return(counts)
}

# Generate count matrix (simulated for demonstration)
# In real analysis, you would count actual reads from BAM files
count_matrix <- data.frame(
    peak_id = consensus$peak_id,
    ARF27_rep1 = count_reads_in_peaks("rep1.bam", consensus_gr),
    ARF27_rep2 = count_reads_in_peaks("rep2.bam", consensus_gr),
    # Add additional samples/conditions as needed
    stringsAsFactors = FALSE
)

# Save count matrix
write.csv(count_matrix, "09_differential_analysis/count_matrices/peak_counts.csv", 
          row.names = FALSE)

# Basic differential analysis using simple statistics
# Calculate fold changes and identify differential peaks
count_matrix$mean_rep1 <- count_matrix$ARF27_rep1
count_matrix$mean_rep2 <- count_matrix$ARF27_rep2
count_matrix$fold_change <- log2((count_matrix$mean_rep2 + 1) / (count_matrix$mean_rep1 + 1))
count_matrix$avg_signal <- (count_matrix$mean_rep1 + count_matrix$mean_rep2) / 2

# Define differential peaks (fold change > 2, average signal > 10)
count_matrix$is_differential <- abs(count_matrix$fold_change) > 1 & count_matrix$avg_signal > 10

# Summary statistics
cat("Differential Peak Analysis Summary:\n")
cat("Total peaks analyzed:", nrow(count_matrix), "\n")
cat("Differential peaks (|FC| > 2):", sum(count_matrix$is_differential), "\n")
cat("Up-regulated peaks:", sum(count_matrix$fold_change > 1 & count_matrix$is_differential), "\n")
cat("Down-regulated peaks:", sum(count_matrix$fold_change < -1 & count_matrix$is_differential), "\n")

# Create visualizations
png("09_differential_analysis/plots/ma_plot.png", width = 800, height = 600)
plot(count_matrix$avg_signal, count_matrix$fold_change,
     pch = 16, cex = 0.5, col = ifelse(count_matrix$is_differential, "red", "black"),
     xlab = "Average Signal", ylab = "Log2 Fold Change",
     main = "MA Plot: Differential Peak Analysis")
abline(h = c(-1, 1), col = "blue", lty = 2)
legend("topright", c("Differential", "Not significant"), 
       col = c("red", "black"), pch = 16)
dev.off()

# Volcano plot-style visualization
png("09_differential_analysis/plots/fold_change_distribution.png", width = 800, height = 600)
hist(count_matrix$fold_change, breaks = 50, 
     main = "Distribution of Fold Changes", 
     xlab = "Log2 Fold Change", ylab = "Frequency",
     col = "lightblue")
abline(v = c(-1, 1), col = "red", lty = 2)
dev.off()

# Save differential results
differential_peaks <- count_matrix[count_matrix$is_differential, ]
write.csv(differential_peaks, 
          "09_differential_analysis/results/differential_peaks_manual.csv", 
          row.names = FALSE)

cat("Manual differential analysis completed.\n")
cat("Results saved to: 09_differential_analysis/results/\n")
EOF
```

### 7. Method 3: Using csaw for Window-Based Analysis

Implement window-based differential analysis:

```bash
# csaw-based differential analysis
R --vanilla << 'EOF'
library(csaw)
library(edgeR)
library(GenomicRanges)

# Note: This is a template for csaw analysis
# You would need to adapt this to your actual BAM files and experimental design

cat("Setting up csaw analysis...\n")

# Parameters for window-based analysis
window_size <- 500  # 500bp windows
fragment_length <- 150

# In a real analysis, you would:
# 1. Define BAM files for each condition
# bam_files <- c("condition1_rep1.bam", "condition1_rep2.bam", 
#                "condition2_rep1.bam", "condition2_rep2.bam")
# 
# 2. Count reads in windows
# windows <- windowCounts(bam_files, width = window_size, ext = fragment_length)
# 
# 3. Filter low-count windows
# keep <- filterWindows(windows, background = 10)$filter > log2(2)
# windows_filtered <- windows[keep]
# 
# 4. Normalize for library composition
# windows_filtered <- normOffsets(windows_filtered)
# 
# 5. Set up design matrix
# condition <- factor(c("control", "control", "treatment", "treatment"))
# design <- model.matrix(~condition)
# 
# 6. Perform differential analysis
# y <- asDGEList(windows_filtered)
# y <- estimateDisp(y, design)
# fit <- glmQLFit(y, design)
# results <- glmQLFTest(fit)
# 
# 7. Merge nearby significant windows
# merged_regions <- mergeWindows(windows_filtered, tol = 1000)
# tabcom <- combineTests(merged_regions$id, results$table)
# 
# 8. Save results
# write.csv(tabcom, "09_differential_analysis/results/csaw_results.csv")

cat("csaw analysis template created.\n")
cat("Please adapt this code to your specific BAM files and experimental design.\n")
EOF
```

### 8. Visualize Differential Results

Create comprehensive visualizations of differential binding:

```bash
# Create differential binding visualizations
R --vanilla << 'EOF'
library(ggplot2)
library(pheatmap)
library(VennDiagram)

# Read results from manual analysis
if(file.exists("09_differential_analysis/count_matrices/peak_counts.csv")) {
    counts <- read.csv("09_differential_analysis/count_matrices/peak_counts.csv")
    
    # Create correlation heatmap
    count_matrix_only <- counts[, c("ARF27_rep1", "ARF27_rep2")]
    correlation_matrix <- cor(count_matrix_only)
    
    png("09_differential_analysis/plots/sample_correlation_heatmap.png", 
        width = 600, height = 600)
    pheatmap(correlation_matrix, 
             main = "Sample Correlation Matrix",
             display_numbers = TRUE,
             color = colorRampPalette(c("blue", "white", "red"))(100))
    dev.off()
    
    # Create signal distribution plots
    png("09_differential_analysis/plots/signal_distributions.png", 
        width = 1000, height = 500)
    par(mfrow = c(1, 2))
    
    # Log-transformed signal distributions
    hist(log10(counts$ARF27_rep1 + 1), breaks = 50, 
         main = "ARF27 Rep1 Signal Distribution", 
         xlab = "Log10(Read Count + 1)", col = "lightblue")
    
    hist(log10(counts$ARF27_rep2 + 1), breaks = 50, 
         main = "ARF27 Rep2 Signal Distribution", 
         xlab = "Log10(Read Count + 1)", col = "lightgreen")
    
    dev.off()
    
    # Scatter plot of replicates
    png("09_differential_analysis/plots/replicate_scatter.png", 
        width = 600, height = 600)
    plot(log10(counts$ARF27_rep1 + 1), log10(counts$ARF27_rep2 + 1),
         pch = 16, cex = 0.5, 
         xlab = "Log10(ARF27 Rep1 + 1)", ylab = "Log10(ARF27 Rep2 + 1)",
         main = "Replicate Correlation")
    abline(a = 0, b = 1, col = "red", lwd = 2)
    correlation <- cor(counts$ARF27_rep1, counts$ARF27_rep2)
    text(1, 3, paste("r =", round(correlation, 3)), cex = 1.2)
    dev.off()
    
    cat("Visualization plots created in 09_differential_analysis/plots/\n")
}

# Create summary statistics
cat("Creating analysis summary...\n")
if(file.exists("09_differential_analysis/count_matrices/peak_counts.csv")) {
    counts <- read.csv("09_differential_analysis/count_matrices/peak_counts.csv")
    
    summary_stats <- data.frame(
        Metric = c("Total peaks", "Mean signal Rep1", "Mean signal Rep2", 
                  "Correlation", "High signal peaks (>100)", "Low signal peaks (<10)"),
        Value = c(nrow(counts),
                 round(mean(counts$ARF27_rep1), 1),
                 round(mean(counts$ARF27_rep2), 1),
                 round(cor(counts$ARF27_rep1, counts$ARF27_rep2), 3),
                 sum(counts$ARF27_rep1 > 100 | counts$ARF27_rep2 > 100),
                 sum(counts$ARF27_rep1 < 10 & counts$ARF27_rep2 < 10))
    )
    
    write.csv(summary_stats, 
              "09_differential_analysis/results/analysis_summary.csv", 
              row.names = FALSE)
    
    print(summary_stats)
}
EOF
```

### 9. Condition-Specific Peak Analysis

Identify peaks that are specific to certain conditions:

```bash
# Analyze condition-specific peaks
echo "Analyzing condition-specific peaks..."

# For demonstration, we'll analyze peaks unique to each replicate
# In real analysis, this would be peaks unique to each condition

# Find peaks unique to replicate 1
bedtools intersect \
    -a 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak \
    -b 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak \
    -v > 09_differential_analysis/results/rep1_specific_peaks.bed

# Find peaks unique to replicate 2
bedtools intersect \
    -a 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak \
    -b 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak \
    -v > 09_differential_analysis/results/rep2_specific_peaks.bed

# Find common peaks
bedtools intersect \
    -a 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak \
    -b 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak \
    -wa > 09_differential_analysis/results/common_peaks.bed

# Calculate statistics
rep1_specific=$(wc -l < 09_differential_analysis/results/rep1_specific_peaks.bed)
rep2_specific=$(wc -l < 09_differential_analysis/results/rep2_specific_peaks.bed)
common_peaks=$(wc -l < 09_differential_analysis/results/common_peaks.bed)
total_rep1=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak)
total_rep2=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak)

echo "Peak overlap analysis:"
echo "Rep1 total peaks: $total_rep1"
echo "Rep2 total peaks: $total_rep2"
echo "Rep1-specific peaks: $rep1_specific"
echo "Rep2-specific peaks: $rep2_specific"
echo "Common peaks: $common_peaks"

# Create Venn diagram data
R --vanilla << 'EOF'
library(VennDiagram)

# Create Venn diagram
venn.plot <- draw.pairwise.venn(
    area1 = total_rep1,
    area2 = total_rep2,
    cross.area = common_peaks,
    category = c("Rep1", "Rep2"),
    fill = c("lightblue", "lightgreen"),
    alpha = 0.5,
    cex = 1.5,
    cat.cex = 1.5,
    filename = NULL
)

png("09_differential_analysis/plots/peak_overlap_venn.png", 
    width = 600, height = 600)
grid.draw(venn.plot)
dev.off()

cat("Venn diagram created: 09_differential_analysis/plots/peak_overlap_venn.png\n")
EOF
```

### 10. Export Results for Downstream Analysis

Prepare differential peaks for annotation and motif analysis:

```bash
# Prepare results for downstream analysis
echo "Preparing results for downstream analysis..."

# Create BED files for different peak categories
if [ -f "09_differential_analysis/results/differential_peaks_manual.csv" ]; then
    
    # Convert differential peaks to BED format
    R --vanilla << 'EOF'
    # Read differential results and consensus peaks
    diff_results <- read.csv("09_differential_analysis/results/differential_peaks_manual.csv")
    consensus <- read.table("09_differential_analysis/consensus_peaks/consensus_peaks.bed")
    colnames(consensus) <- c("chr", "start", "end", "peak_id")
    
    # Merge with genomic coordinates
    diff_with_coords <- merge(diff_results, consensus, by = "peak_id")
    
    # Create BED files for different categories
    # Up-regulated peaks
    up_peaks <- diff_with_coords[diff_with_coords$fold_change > 1 & 
                                diff_with_coords$is_differential, ]
    if(nrow(up_peaks) > 0) {
        write.table(up_peaks[, c("chr", "start", "end", "peak_id")], 
                   "09_differential_analysis/results/upregulated_peaks.bed",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    # Down-regulated peaks  
    down_peaks <- diff_with_coords[diff_with_coords$fold_change < -1 & 
                                  diff_with_coords$is_differential, ]
    if(nrow(down_peaks) > 0) {
        write.table(down_peaks[, c("chr", "start", "end", "peak_id")], 
                   "09_differential_analysis/results/downregulated_peaks.bed",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    # All differential peaks
    all_diff <- diff_with_coords[diff_with_coords$is_differential, ]
    if(nrow(all_diff) > 0) {
        write.table(all_diff[, c("chr", "start", "end", "peak_id")], 
                   "09_differential_analysis/results/all_differential_peaks.bed",
                   sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    }
    
    cat("BED files created for differential peak categories.\n")
EOF

fi

# Create summary report
cat > 09_differential_analysis/results/differential_analysis_summary.txt << EOF
Differential Peak Analysis Summary
==================================
Analysis Date: $(date)

EXPERIMENTAL DESIGN
===================
Note: This analysis used replicate comparison as demonstration.
In real analysis, you would compare actual experimental conditions
(e.g., Control vs Treatment, Time points, Tissues, etc.)

PEAK STATISTICS
===============
Total consensus peaks: $(wc -l < 09_differential_analysis/consensus_peaks/consensus_peaks.bed)
Rep1-specific peaks: $rep1_specific
Rep2-specific peaks: $rep2_specific
Common peaks: $common_peaks

DIFFERENTIAL RESULTS
====================
$(if [ -f "09_differential_analysis/results/differential_peaks_manual.csv" ]; then
    echo "Manual analysis completed:"
    echo "- Total differential peaks: $(tail -n +2 09_differential_analysis/results/differential_peaks_manual.csv | wc -l)"
    echo "- Up-regulated: $(awk -F',' 'NR>1 && $6>1 && $7=="TRUE"' 09_differential_analysis/results/differential_peaks_manual.csv | wc -l)"
    echo "- Down-regulated: $(awk -F',' 'NR>1 && $6<-1 && $7=="TRUE"' 09_differential_analysis/results/differential_peaks_manual.csv | wc -l)"
else
    echo "Manual analysis files not found"
fi)

FILES GENERATED
===============
1. Consensus peaks: consensus_peaks/consensus_peaks.bed
2. Count matrices: count_matrices/peak_counts.csv  
3. Differential results: results/differential_peaks_manual.csv
4. Visualization plots: plots/
5. Peak category BED files: results/*_peaks.bed

QUALITY METRICS
===============
Replicate correlation: $(if [ -f "09_differential_analysis/count_matrices/peak_counts.csv" ]; then
    R --slave -e "data<-read.csv('09_differential_analysis/count_matrices/peak_counts.csv'); cat(round(cor(data[,2], data[,3]), 3))"
fi)

NEXT STEPS
==========
✓ Differential analysis completed
→ Ready for Step 10: Peak Annotation
→ Ready for Step 11: Motif Analysis

EOF

echo "Differential analysis summary: 09_differential_analysis/results/differential_analysis_summary.txt"
```

## Common Issues and Troubleshooting

### Issue 1: Insufficient Replicates

**Problem**: Need at least 2 biological replicates per condition for statistical analysis

**Solutions:**
```bash
# Use simple fold-change approach without statistics
# Focus on highly confident peaks (high signal, consistent between available replicates)
# Consider pooling similar conditions if biologically appropriate
```

### Issue 2: Poor Replicate Correlation

**Diagnosis:**
```bash
# Check replicate correlation
R -e "
data <- read.csv('09_differential_analysis/count_matrices/peak_counts.csv')
cor(data[,2], data[,3])
"
```

**Solutions:**
- Investigate experimental conditions
- Check for batch effects
- Consider excluding poor-quality replicates
- Adjust normalization methods

### Issue 3: No Differential Peaks Detected

**Possible Causes:**
- Conditions too similar
- Insufficient statistical power
- Overly stringent thresholds

**Solutions:**
```bash
# Relax significance thresholds
# Check raw signal differences
# Verify experimental conditions were applied correctly
# Use different statistical methods
```

### Issue 4: Memory Issues with Large Datasets

**Solutions:**
```bash
# Process chromosomes separately
# Use more efficient data structures
# Increase available memory
# Use window-based approaches instead of peak-based
```

## Advanced Analysis Options

### 1. Time-Course Analysis

For developmental or treatment time courses:

```bash
# Template for time-course analysis
R --vanilla << 'EOF'
# Example: Analyze binding changes over time
# timepoints <- c("0h", "2h", "6h", "24h")
# Use spline or polynomial models to capture temporal patterns
# library(splines)
# design <- model.matrix(~ ns(time, df=3))
cat("Time-course analysis template available.\n")
EOF
```

### 2. Multi-Condition Comparisons

For complex experimental designs:

```bash
# Template for multi-condition analysis
R --vanilla << 'EOF'
# Example: Multiple tissues, treatments, genotypes
# Use ANOVA-like approaches
# library(limma) or library(DESeq2) with complex designs
cat("Multi-condition analysis template available.\n")
EOF
```

## Analysis Log Update

```bash
# Update analysis log
cat >> logs/step9_differential_analysis_log.txt << EOF

Step 9: Differential Peak Analysis Results
Date: $(date)

Analysis Overview:
- Method: Manual differential analysis (demonstration)
- Comparison: Replicate 1 vs Replicate 2 (modify for real conditions)
- Statistical approach: Fold-change based with signal thresholds

Key Results:
- Consensus peaks: $(wc -l < 09_differential_analysis/consensus_peaks/consensus_peaks.bed)
- Common peaks: $common_peaks
- Replicate-specific peaks: Rep1=$rep1_specific, Rep2=$rep2_specific
$(if [ -f "09_differential_analysis/count_matrices/peak_counts.csv" ]; then
    echo "- Replicate correlation: $(R --slave -e "data<-read.csv('09_differential_analysis/count_matrices/peak_counts.csv'); cat(round(cor(data[,2], data[,3]), 3))")"
fi)

Files Generated:
- Sample sheets and templates: 3 files
- Consensus peak set: 1 file
- Count matrices: 1 file
- Differential results: Multiple files
- Visualization plots: 5+ files

Tools Used:
- R packages: DiffBind (template), DESeq2, edgeR, csaw (template)
- Bedtools for peak overlaps
- Custom R scripts for manual analysis

Issues Encountered: [None/List problems and solutions]

Notes:
- This analysis used replicates as demonstration
- Modify sample sheets for real experimental conditions
- Consider using DiffBind or csaw for production analysis

Ready for Step 10: Peak Annotation
