# Step 8: Peak Data Quality Control for ChIP-seq

## Objective
Perform comprehensive quality assessment of filtered ARF27 peaks using established ChIP-seq metrics to ensure high-quality binding site predictions and validate experimental success before downstream analysis.

## Why Peak Quality Control is Critical

### Quality Control Validates:
- **Experimental success**: ChIP efficiency and specificity
- **Peak reliability**: Signal-to-noise ratios and statistical significance
- **Biological relevance**: Expected binding patterns and genomic distribution
- **Technical quality**: Consistency between replicates and appropriate controls

### ChIP-seq Quality Metrics:
1. **FRiP Score**: Fraction of Reads in Peaks
2. **NSC/RSC**: Normalized/Relative Strand Cross-correlation
3. **Peak enrichment**: Signal intensity and fold-change
4. **Fragment length distribution**: Insert size patterns
5. **Peak annotation**: Genomic feature distribution
6. **Motif enrichment**: Expected transcription factor motifs

### ARF27-Specific Expectations:
- **FRiP score**: 5-30% for transcription factors
- **Peak location**: Promoters and enhancers of auxin-responsive genes
- **Motif presence**: Auxin Response Elements (AREs)
- **Peak sharpness**: Narrow peaks with defined summits

## Step-by-Step Quality Control Analysis

### 1. Set Up QC Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create QC directory structure
mkdir -p 08_peak_qc
mkdir -p 08_peak_qc/frip_analysis
mkdir -p 08_peak_qc/fragment_analysis
mkdir -p 08_peak_qc/enrichment_analysis
mkdir -p 08_peak_qc/cross_correlation
mkdir -p 08_peak_qc/genomic_distribution
mkdir -p 08_peak_qc/summary_plots
```

### 2. Calculate FRiP Scores (Fraction of Reads in Peaks)

FRiP is a key metric indicating ChIP enrichment quality:

```bash
# Calculate FRiP scores for all samples
echo "Calculating FRiP scores..."
echo "Sample\tTotal_Reads\tReads_in_Peaks\tFRiP_Score" > 08_peak_qc/frip_analysis/frip_scores.txt

# Function to calculate FRiP
calculate_frip() {
    local bam_file=$1
    local peak_file=$2
    local sample_name=$3
    
    # Count total reads
    total_reads=$(samtools view -c $bam_file)
    
    # Count reads overlapping peaks
    reads_in_peaks=$(bedtools intersect -a $bam_file -b $peak_file -u | samtools view -c -)
    
    # Calculate FRiP score
    frip_score=$(echo "scale=4; $reads_in_peaks * 100 / $total_reads" | bc)
    
    echo -e "${sample_name}\t${total_reads}\t${reads_in_peaks}\t${frip_score}" >> 08_peak_qc/frip_analysis/frip_scores.txt
    echo "$sample_name FRiP: $frip_score%"
}

# Calculate FRiP for individual replicates
calculate_frip "05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam" \
               "07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak" \
               "ARF27_ChIP_rep1"

calculate_frip "05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam" \
               "07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak" \
               "ARF27_ChIP_rep2"

# Calculate FRiP for merged data
calculate_frip "06_peaks/merged_replicates/ARF27_ChIP_merged.bam" \
               "07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak" \
               "ARF27_ChIP_merged"

# Calculate FRiP for input controls (should be low)
calculate_frip "05_processed_bams/deduplicated/Input_control_rep1_dedup.bam" \
               "07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak" \
               "Input_control_rep1"

calculate_frip "05_processed_bams/deduplicated/Input_control_rep2_dedup.bam" \
               "07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak" \
               "Input_control_rep2"

# Display FRiP summary
echo "FRiP Score Summary:"
column -t 08_peak_qc/frip_analysis/frip_scores.txt
```

### 3. Fragment Length Analysis

Analyze fragment size distribution to assess ChIP quality:

```bash
# Calculate insert size statistics
echo "Analyzing fragment length distribution..."

# Extract properly paired reads and calculate insert sizes
samtools view -f 2 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam | \
awk '$9 > 0 {print $9}' > 08_peak_qc/fragment_analysis/ARF27_rep1_insert_sizes.txt

samtools view -f 2 05_processed_bams/deduplicated/ARF27_ChIP_rep2_dedup.bam | \
awk '$9 > 0 {print $9}' > 08_peak_qc/fragment_analysis/ARF27_rep2_insert_sizes.txt

# Analyze fragment length statistics
R --vanilla << 'EOF'
# Read fragment length data
rep1_sizes <- read.table("08_peak_qc/fragment_analysis/ARF27_rep1_insert_sizes.txt")$V1
rep2_sizes <- read.table("08_peak_qc/fragment_analysis/ARF27_rep2_insert_sizes.txt")$V1

# Filter realistic fragment sizes (50-1000 bp)
rep1_filtered <- rep1_sizes[rep1_sizes >= 50 & rep1_sizes <= 1000]
rep2_filtered <- rep2_sizes[rep2_sizes >= 50 & rep2_sizes <= 1000]

# Calculate statistics
cat("Fragment Length Statistics:\n")
cat("Rep1 - Mean:", round(mean(rep1_filtered)), "bp, Median:", round(median(rep1_filtered)), "bp\n")
cat("Rep2 - Mean:", round(mean(rep2_filtered)), "bp, Median:", round(median(rep2_filtered)), "bp\n")

# Create fragment length distribution plot
png("08_peak_qc/fragment_analysis/fragment_length_distribution.png", width=1000, height=600)
par(mfrow=c(1,2))

hist(rep1_filtered, breaks=50, main="ARF27 Rep1 Fragment Lengths", 
     xlab="Fragment Length (bp)", ylab="Frequency", col="lightblue", xlim=c(50, 500))
abline(v=median(rep1_filtered), col="red", lwd=2, lty=2)

hist(rep2_filtered, breaks=50, main="ARF27 Rep2 Fragment Lengths", 
     xlab="Fragment Length (bp)", ylab="Frequency", col="lightgreen", xlim=c(50, 500))
abline(v=median(rep2_filtered), col="red", lwd=2, lty=2)

dev.off()

# Save statistics
sink("08_peak_qc/fragment_analysis/fragment_stats.txt")
cat("Fragment Length Analysis\n")
cat("========================\n\n")
cat("ARF27 ChIP Rep1:\n")
cat("  Count:", length(rep1_filtered), "fragments\n")
cat("  Mean:", round(mean(rep1_filtered)), "bp\n")
cat("  Median:", round(median(rep1_filtered)), "bp\n")
cat("  SD:", round(sd(rep1_filtered)), "bp\n")
cat("  Range:", min(rep1_filtered), "-", max(rep1_filtered), "bp\n\n")

cat("ARF27 ChIP Rep2:\n")
cat("  Count:", length(rep2_filtered), "fragments\n")
cat("  Mean:", round(mean(rep2_filtered)), "bp\n")
cat("  Median:", round(median(rep2_filtered)), "bp\n")
cat("  SD:", round(sd(rep2_filtered)), "bp\n")
cat("  Range:", min(rep2_filtered), "-", max(rep2_filtered), "bp\n")
sink()

# Expected: Mean ~150-200bp for good ChIP-seq
EOF
```

### 4. Cross-Correlation Analysis

Calculate NSC and RSC metrics using phantompeakqualtools:

```bash
# Install phantompeakqualtools if not available
# R -e "install.packages(c('caTools', 'snow', 'snowfall', 'bitops', 'RCurl'))"
# R -e "source('http://phantompeakqualtools.googlecode.com/files/phantompeakqualtools.R')"

# Alternative: Use deepTools for cross-correlation analysis
echo "Calculating cross-correlation metrics..."

# Calculate cross-correlation for merged ChIP sample
plotFingerprint \
    --bamfiles 06_peaks/merged_replicates/ARF27_ChIP_merged.bam \
               06_peaks/merged_replicates/Input_control_merged.bam \
    --labels "ARF27_ChIP" "Input_Control" \
    --plotFile 08_peak_qc/cross_correlation/fingerprint_plot.png \
    --outRawCounts 08_peak_qc/cross_correlation/fingerprint_data.txt \
    --numberOfProcessors 8 \
    --plotTitle "ARF27 ChIP-seq Fingerprint Plot"

# Calculate enrichment at TSS
computeMatrix reference-point \
    --referencePoint TSS \
    --beforeRegionStartLength 2000 \
    --afterRegionStartLength 2000 \
    --scoreFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_merged_RPM.bw \
                   05_processed_bams/normalized_bigwigs/Input_control_merged_RPM.bw \
    --regionsFileName ../reference_genome/maize_genes_TSS.bed \
    --out 08_peak_qc/cross_correlation/TSS_matrix.gz \
    --numberOfProcessors 8

# Create TSS enrichment plot
plotHeatmap \
    --matrixFile 08_peak_qc/cross_correlation/TSS_matrix.gz \
    --plotFile 08_peak_qc/cross_correlation/TSS_enrichment.png \
    --plotTitle "ARF27 Enrichment at TSS"
```

### 5. Peak Signal Enrichment Analysis

Analyze peak signal strength and fold-enrichment:

```bash
# Extract peak signal statistics
echo "Analyzing peak signal enrichment..."

# Extract signal values from narrowPeak files
awk '{print $7}' 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak > \
    08_peak_qc/enrichment_analysis/peak_signals.txt

awk '{print $8}' 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak > \
    08_peak_qc/enrichment_analysis/peak_pvalues.txt

awk '{print $9}' 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak > \
    08_peak_qc/enrichment_analysis/peak_qvalues.txt

# Analyze signal distribution
R --vanilla << 'EOF'
# Read peak statistics
signals <- read.table("08_peak_qc/enrichment_analysis/peak_signals.txt")$V1
pvalues <- read.table("08_peak_qc/enrichment_analysis/peak_pvalues.txt")$V1
qvalues <- read.table("08_peak_qc/enrichment_analysis/peak_qvalues.txt")$V1

# Calculate enrichment statistics
cat("Peak Enrichment Statistics:\n")
cat("Signal fold-change:\n")
cat("  Mean:", round(mean(signals), 2), "\n")
cat("  Median:", round(median(signals), 2), "\n")
cat("  >2-fold peaks:", sum(signals > 2), "/", length(signals), 
    "(", round(100*sum(signals > 2)/length(signals), 1), "%)\n")
cat("  >5-fold peaks:", sum(signals > 5), "/", length(signals), 
    "(", round(100*sum(signals > 5)/length(signals), 1), "%)\n")

cat("\nP-value distribution:\n")
cat("  Mean -log10(p):", round(mean(pvalues), 2), "\n")
cat("  Median -log10(p):", round(median(pvalues), 2), "\n")

# Create enrichment plots
png("08_peak_qc/enrichment_analysis/peak_enrichment_plots.png", width=1200, height=800)
par(mfrow=c(2,2))

# Signal fold-change distribution
hist(signals, breaks=50, main="Peak Signal Distribution", 
     xlab="Fold Enrichment", ylab="Frequency", col="lightblue")
abline(v=median(signals), col="red", lwd=2, lty=2)

# Log-transformed signal for better visualization
hist(log2(signals), breaks=50, main="Peak Signal Distribution (log2)", 
     xlab="Log2 Fold Enrichment", ylab="Frequency", col="lightgreen")
abline(v=log2(median(signals)), col="red", lwd=2, lty=2)

# P-value distribution
hist(pvalues, breaks=50, main="Peak P-value Distribution", 
     xlab="-log10(p-value)", ylab="Frequency", col="lightcoral")

# Q-value distribution
hist(qvalues, breaks=50, main="Peak Q-value Distribution", 
     xlab="-log10(q-value)", ylab="Frequency", col="lightyellow")

dev.off()

# Save summary
sink("08_peak_qc/enrichment_analysis/enrichment_summary.txt")
cat("Peak Signal Enrichment Summary\n")
cat("==============================\n\n")
cat("Total peaks analyzed:", length(signals), "\n")
cat("Signal fold-change statistics:\n")
cat("  Mean:", round(mean(signals), 2), "\n")
cat("  Median:", round(median(signals), 2), "\n")
cat("  Standard deviation:", round(sd(signals), 2), "\n")
cat("  Range:", round(min(signals), 2), "-", round(max(signals), 2), "\n")
cat("\nEnrichment categories:\n")
cat("  >2-fold:", sum(signals > 2), "(", round(100*sum(signals > 2)/length(signals), 1), "%)\n")
cat("  >5-fold:", sum(signals > 5), "(", round(100*sum(signals > 5)/length(signals), 1), "%)\n")
cat("  >10-fold:", sum(signals > 10), "(", round(100*sum(signals > 10)/length(signals), 1), "%)\n")
sink()
EOF
```

### 6. Genomic Distribution Analysis

Analyze where peaks are located relative to genes and genomic features:

```bash
# Download gene annotation if not already available
# (Assuming you have a GFF3 file for maize B73-v4)

# Convert GFF3 to BED for different features
echo "Analyzing genomic distribution of peaks..."

# Extract different genomic features (if annotation available)
if [ -f "../reference_genome/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3" ]; then
    
    # Extract promoter regions (2kb upstream of TSS)
    awk '$3=="gene"' ../reference_genome/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 | \
    awk '{if($7=="+") print $1"\t"($4-2000)"\t"$4"\tpromoter"; 
          else print $1"\t"$5"\t"($5+2000)"\tpromoter"}' | \
    awk '$2>=0' > 08_peak_qc/genomic_distribution/promoter_regions.bed
    
    # Extract gene bodies
    awk '$3=="gene" {print $1"\t"($4-1)"\t"$5"\tgene_body"}' \
        ../reference_genome/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 > \
        08_peak_qc/genomic_distribution/gene_bodies.bed
    
    # Extract intergenic regions (simplified approach)
    # This is a placeholder - proper intergenic requires more complex processing
    
    echo "Calculating peak overlap with genomic features..."
    
    # Count peaks in different features
    promoter_peaks=$(bedtools intersect \
        -a 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak \
        -b 08_peak_qc/genomic_distribution/promoter_regions.bed \
        -wa | wc -l)
    
    gene_body_peaks=$(bedtools intersect \
        -a 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak \
        -b 08_peak_qc/genomic_distribution/gene_bodies.bed \
        -wa | wc -l)
    
    total_peaks=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
    
    # Calculate percentages
    promoter_percent=$(echo "scale=1; $promoter_peaks * 100 / $total_peaks" | bc)
    gene_body_percent=$(echo "scale=1; $gene_body_peaks * 100 / $total_peaks" | bc)
    other_percent=$(echo "scale=1; 100 - $promoter_percent - $gene_body_percent" | bc)
    
    echo "Genomic Distribution of ARF27 Peaks:"
    echo "Promoters: $promoter_peaks ($promoter_percent%)"
    echo "Gene bodies: $gene_body_peaks ($gene_body_percent%)"
    echo "Other/Intergenic: $other_percent%"
    
    # Create distribution summary
    cat > 08_peak_qc/genomic_distribution/distribution_summary.txt << EOF
Genomic Distribution of ARF27 Peaks
====================================

Total peaks: $total_peaks

Distribution:
  Promoters (2kb upstream): $promoter_peaks ($promoter_percent%)
  Gene bodies: $gene_body_peaks ($gene_body_percent%)
  Other/Intergenic: $other_percent%

Expected for transcription factors:
  Promoters: 30-60%
  Gene bodies: 20-40%
  Intergenic: 20-40%
EOF

else
    echo "Gene annotation file not found. Skipping genomic distribution analysis."
    echo "Please provide: ../reference_genome/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3"
fi
```

### 7. Peak Width and Summit Analysis

Analyze peak characteristics in detail:

```bash
# Analyze peak widths and summit positions
echo "Analyzing peak characteristics..."

# Extract peak information
awk '{print $3-$2, $10, $5}' 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak > \
    08_peak_qc/summary_plots/peak_characteristics.txt

R --vanilla << 'EOF'
# Read peak characteristics
data <- read.table("08_peak_qc/summary_plots/peak_characteristics.txt")
colnames(data) <- c("width", "summit_offset", "score")

# Calculate statistics
cat("Peak Characteristics Analysis:\n")
cat("Peak widths:\n")
cat("  Mean:", round(mean(data$width)), "bp\n")
cat("  Median:", round(median(data$width)), "bp\n")
cat("  Standard deviation:", round(sd(data$width)), "bp\n")

cat("\nSummit positions (relative to peak start):\n")
cat("  Mean offset:", round(mean(data$summit_offset)), "bp\n")
cat("  Median offset:", round(median(data$summit_offset)), "bp\n")

# Create comprehensive peak characteristics plot
png("08_peak_qc/summary_plots/peak_characteristics.png", width=1200, height=800)
par(mfrow=c(2,2))

# Peak width distribution
hist(data$width, breaks=50, main="Peak Width Distribution", 
     xlab="Peak Width (bp)", ylab="Frequency", col="lightblue")
abline(v=median(data$width), col="red", lwd=2, lty=2)

# Summit position distribution
hist(data$summit_offset, breaks=50, main="Summit Position Distribution", 
     xlab="Summit Offset from Peak Start (bp)", ylab="Frequency", col="lightgreen")
abline(v=median(data$summit_offset), col="red", lwd=2, lty=2)

# Score vs width relationship
plot(data$width, data$score, pch=16, cex=0.5, 
     main="Peak Score vs Width", xlab="Peak Width (bp)", ylab="MACS2 Score")

# Score distribution
hist(data$score, breaks=50, main="Peak Score Distribution", 
     xlab="MACS2 Score", ylab="Frequency", col="lightcoral")

dev.off()

# Save detailed statistics
sink("08_peak_qc/summary_plots/peak_characteristics_summary.txt")
cat("Peak Characteristics Summary\n")
cat("============================\n\n")
cat("Peak Width Statistics:\n")
cat("  Count:", nrow(data), "peaks\n")
cat("  Mean:", round(mean(data$width)), "bp\n")
cat("  Median:", round(median(data$width)), "bp\n")
cat("  Standard deviation:", round(sd(data$width)), "bp\n")
cat("  Range:", min(data$width), "-", max(data$width), "bp\n")
cat("  Quartiles:", quantile(data$width), "\n\n")

cat("Summit Position Statistics:\n")
cat("  Mean offset:", round(mean(data$summit_offset)), "bp\n")
cat("  Median offset:", round(median(data$summit_offset)), "bp\n")
cat("  Standard deviation:", round(sd(data$summit_offset)), "bp\n")

cat("\nPeak Quality Categories:\n")
cat("  Narrow peaks (<200bp):", sum(data$width < 200), 
    "(", round(100*sum(data$width < 200)/nrow(data), 1), "%)\n")
cat("  Medium peaks (200-500bp):", sum(data$width >= 200 & data$width <= 500), 
    "(", round(100*sum(data$width >= 200 & data$width <= 500)/nrow(data), 1), "%)\n")
cat("  Broad peaks (>500bp):", sum(data$width > 500), 
    "(", round(100*sum(data$width > 500)/nrow(data), 1), "%)\n")
sink()
EOF
```

### 8. Peak Shape Analysis

```bash
# Analyze peak shapes using deepTools
echo "Analyzing peak shapes and profiles..."

# Create matrix for peak shape analysis
computeMatrix reference-point \
    --referencePoint center \
    --beforeRegionStartLength 500 \
    --afterRegionStartLength 500 \
    --scoreFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_merged_RPM.bw \
                   05_processed_bams/normalized_bigwigs/Input_control_merged_RPM.bw \
    --regionsFileName 07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed \
    --out 08_peak_qc/peak_shape_matrix.gz \
    --numberOfProcessors 8

# Generate peak shape heatmap
plotHeatmap \
    --matrixFile 08_peak_qc/peak_shape_matrix.gz \
    --plotFile 08_peak_qc/peak_shape_heatmap.png \
    --plotTitle "ARF27 Peak Shape Analysis" \
    --colorMap Blues \
    --sortRegions descend

# Generate average profile
plotProfile \
    --matrixFile 08_peak_qc/peak_shape_matrix.gz \
    --plotFile 08_peak_qc/peak_average_profile.png \
    --plotTitle "Average ARF27 Peak Profile" \
    --colors blue red \
    --samplesLabel "ARF27_ChIP" "Input_Control"
```

### 9. Peak Reproducibility Analysis (IDR)

```bash
# Install IDR if not available
# pip install idr

# Prepare peaks for IDR analysis
echo "Preparing peaks for IDR analysis..."

# Convert narrowPeak to IDR format (only if IDR is available)
# Sort peaks by p-value (column 8)
sort -k8,8nr 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak > \
    08_peak_qc/ARF27_rep1_for_IDR.narrowPeak

sort -k8,8nr 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak > \
    08_peak_qc/ARF27_rep2_for_IDR.narrowPeak

# Calculate basic overlap statistics as IDR alternative
echo "Calculating replicate peak overlap (IDR alternative)..."

# Find peaks present in both replicates
bedtools intersect \
    -a 08_peak_qc/ARF27_rep1_for_IDR.narrowPeak \
    -b 08_peak_qc/ARF27_rep2_for_IDR.narrowPeak \
    -wo > 08_peak_qc/overlapping_peaks.txt

# Calculate overlap statistics
total_rep1=$(wc -l < 08_peak_qc/ARF27_rep1_for_IDR.narrowPeak)
total_rep2=$(wc -l < 08_peak_qc/ARF27_rep2_for_IDR.narrowPeak)
overlapping=$(wc -l < 08_peak_qc/overlapping_peaks.txt)

overlap_rate_rep1=$(echo "scale=1; $overlapping * 100 / $total_rep1" | bc)
overlap_rate_rep2=$(echo "scale=1; $overlapping * 100 / $total_rep2" | bc)

echo "Peak Reproducibility Analysis:"
echo "Rep1 peaks: $total_rep1"
echo "Rep2 peaks: $total_rep2"
echo "Overlapping peaks: $overlapping"
echo "Rep1 overlap rate: $overlap_rate_rep1%"
echo "Rep2 overlap rate: $overlap_rate_rep2%"

# Save reproducibility summary
cat > 08_peak_qc/reproducibility_summary.txt << EOF
Peak Reproducibility Analysis
=============================

Total peaks:
- Replicate 1: $total_rep1
- Replicate 2: $total_rep2

Overlap analysis:
- Overlapping peaks: $overlapping
- Rep1 overlap rate: $overlap_rate_rep1%
- Rep2 overlap rate: $overlap_rate_rep2%

Quality assessment:
- Excellent: >80% overlap
- Good: 60-80% overlap
- Acceptable: 40-60% overlap
- Poor: <40% overlap

Current status: $(if (( $(echo "$overlap_rate_rep1 >= 80" | bc -l) )); then echo "Excellent"; elif (( $(echo "$overlap_rate_rep1 >= 60" | bc -l) )); then echo "Good"; elif (( $(echo "$overlap_rate_rep1 >= 40" | bc -l) )); then echo "Acceptable"; else echo "Poor"; fi)
EOF
```

### 10. Generate Comprehensive QC Report

Create a summary report with all quality metrics:

```bash
# Generate comprehensive QC report
echo "Generating comprehensive QC report..."

# Read various statistics files
frip_data=$(tail -n +2 08_peak_qc/frip_analysis/frip_scores.txt)
fragment_stats=$(cat 08_peak_qc/fragment_analysis/fragment_stats.txt 2>/dev/null || echo "Fragment analysis not available")
enrichment_stats=$(cat 08_peak_qc/enrichment_analysis/enrichment_summary.txt 2>/dev/null || echo "Enrichment analysis not available")
distribution_stats=$(cat 08_peak_qc/genomic_distribution/distribution_summary.txt 2>/dev/null || echo "Distribution analysis not available")
peak_stats=$(cat 08_peak_qc/summary_plots/peak_characteristics_summary.txt 2>/dev/null || echo "Peak characteristics not available")

# Create comprehensive report
cat > 08_peak_qc/ARF27_ChIP_seq_QC_report.txt << EOF
ARF27 ChIP-seq Quality Control Report
=====================================
Date: $(date)
Analysis: Maize B73-v4 genome

SUMMARY METRICS
===============

Peak Counts:
- ARF27 Rep1: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak) peaks
- ARF27 Rep2: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak) peaks
- Merged: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak) peaks

FRiP Scores (Fraction of Reads in Peaks):
$frip_data

DETAILED ANALYSIS
================

$fragment_stats

$enrichment_stats

$distribution_stats

$peak_stats

QUALITY ASSESSMENT
==================

FRiP Score Interpretation:
- ChIP samples: Target >5% for TFs, >10% excellent
- Input controls: Should be <2%

Fragment Length Interpretation:
- Expected: 147bp (nucleosome core) ± 50bp
- Actual: See fragment analysis above

Peak Characteristics Interpretation:
- TF peaks: Typically 100-500bp wide
- Sharp summits: Indicate good specificity
- Score distribution: Higher scores = more confident peaks

RECOMMENDATIONS
===============

Based on the quality metrics:
1. FRiP scores: [Interpret based on actual values]
2. Fragment lengths: [Interpret based on actual values]
3. Peak characteristics: [Interpret based on actual values]
4. Overall assessment: [Good/Acceptable/Poor]

Next steps: [Recommendations based on quality]

EOF

echo "QC report generated: 08_peak_qc/ARF27_ChIP_seq_QC_report.txt"
```

### 11. Quality Control Checklist

Create a final checklist to ensure all quality metrics are acceptable:

```bash
# Generate automated quality assessment
echo "Performing automated quality assessment..."

# Read key metrics
frip_chip=$(grep "ARF27_ChIP_merged" 08_peak_qc/frip_analysis/frip_scores.txt | cut -f4)
total_peaks=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)

# Create quality checklist
cat > 08_peak_qc/quality_checklist.txt << EOF
ARF27 ChIP-seq Quality Control Checklist
=========================================

□ 1. Peak Count Assessment
   Total peaks: $total_peaks
   Expected range: 500-5000 for TF ChIP-seq
   Status: $(if [ "$total_peaks" -ge 500 ] && [ "$total_peaks" -le 5000 ]; then echo "✓ PASS"; else echo "⚠ REVIEW"; fi)

□ 2. FRiP Score Assessment  
   ChIP FRiP: $frip_chip%
   Expected: >5% for TFs
   Status: $(if (( $(echo "$frip_chip >= 5" | bc -l) )); then echo "✓ PASS"; else echo "⚠ REVIEW"; fi)

□ 3. Replicate Consistency
   Overlap rate: $overlap_rate_rep1%
   Expected: >50% for TFs
   Status: $(if (( $(echo "$overlap_rate_rep1 >= 50" | bc -l) )); then echo "✓ PASS"; else echo "⚠ REVIEW"; fi)

□ 4. Fragment Length Distribution
   Check: fragment_analysis/fragment_stats.txt
   Expected: ~147bp ± 50bp
   Status: Manual review required

□ 5. Peak Signal Enrichment
   Check: enrichment_analysis/enrichment_summary.txt
   Expected: >2-fold enrichment for majority of peaks
   Status: Manual review required

□ 6. Genomic Distribution
   Check: genomic_distribution/distribution_summary.txt
   Expected: Promoter enrichment for TFs
   Status: Manual review required

OVERALL ASSESSMENT
==================
Automated checks: $(if (( $(echo "$frip_chip >= 5" | bc -l) )) && [ "$total_peaks" -ge 500 ] && [ "$total_peaks" -le 5000 ] && (( $(echo "$overlap_rate_rep1 >= 50" | bc -l) )); then echo "PASS - Proceed to next step"; else echo "REVIEW REQUIRED - Check failed metrics"; fi)

Manual review items:
- Fragment length distribution
- Peak signal enrichment patterns  
- Genomic distribution patterns
- Visual inspection of browser tracks

Recommendation: 
$(if (( $(echo "$frip_chip >= 5" | bc -l) )) && [ "$total_peaks" -ge 500 ] && [ "$total_peaks" -le 5000 ] && (( $(echo "$overlap_rate_rep1 >= 50" | bc -l) )); then echo "Quality metrics look good. Proceed to differential analysis."; else echo "Some quality metrics need attention. Review failed items before proceeding."; fi)
EOF

echo "Quality checklist created: 08_peak_qc/quality_checklist.txt"
```

## Quality Interpretation Guidelines

### 1. FRiP Score Interpretation

```bash
# Interpret FRiP scores automatically
R --vanilla << 'EOF'
# Read FRiP data
frip_data <- read.table("08_peak_qc/frip_analysis/frip_scores.txt", header=TRUE, sep="\t")

# Define quality thresholds for transcription factors
interpret_frip <- function(sample_name, frip_score) {
    if(grepl("ChIP", sample_name)) {
        if(frip_score >= 15) return("Excellent")
        else if(frip_score >= 10) return("Good")
        else if(frip_score >= 5) return("Acceptable")
        else return("Poor")
    } else if(grepl("Input", sample_name)) {
        if(frip_score <= 2) return("Good")
        else if(frip_score <= 5) return("Acceptable")
        else return("Poor - possible contamination")
    }
}

cat("FRiP Score Quality Assessment:\n")
for(i in 1:nrow(frip_data)) {
    assessment <- interpret_frip(frip_data$Sample[i], frip_data$FRiP_Score[i])
    cat(frip_data$Sample[i], ":", frip_data$FRiP_Score[i], "% -", assessment, "\n")
}
EOF
```

### 2. Expected Quality Benchmarks

| Metric | Excellent | Good | Acceptable | Poor |
|--------|-----------|------|------------|------|
| FRiP (ChIP) | >15% | 10-15% | 5-10% | <5% |
| FRiP (Input) | <1% | 1-2% | 2-5% | >5% |
| Fragment Length | 130-170bp | 100-200bp | 80-250bp | Outside range |
| Peak Width | 100-300bp | 100-500bp | 50-800bp | >1000bp |
| Fold Enrichment | >5x | 3-5x | 2-3x | <2x |

## Common Issues and Troubleshooting

### Issue 1: Low FRiP Scores (<5%)

**Possible Causes:**
- Poor ChIP efficiency
- Inappropriate peak calling parameters
- Low sequencing depth
- Poor antibody quality

**Diagnosis:**
```bash
# Check peak calling parameters
grep "INFO" logs/ARF27_ChIP_merged_peak_calling.log

# Check sequencing depth
samtools flagstat 05_processed_bams/deduplicated/ARF27_ChIP_rep1_dedup.bam
```

**Solutions:**
- Relax peak calling parameters (higher q-value threshold)
- Increase sequencing depth
- Optimize ChIP protocol
- Validate antibody specificity

### Issue 2: High Input FRiP Scores (>5%)

**Possible Causes:**
- Cross-contamination
- Poor input control
- Non-specific enrichment

**Solutions:**
```bash
# Check for contamination patterns
bedtools intersect -a ChIP_peaks.narrowPeak -b Input_peaks.narrowPeak -wa | wc -l

# Use stricter peak calling for input
macs2 callpeak --treatment Input.bam --qvalue 0.001 --name Input_strict
```

### Issue 3: Unusual Fragment Length Distribution

**Expected**: Single peak around 147bp

**Problems:**
- Multiple peaks: Suggest nucleosome laddering
- Very broad distribution: Poor sonication
- Shifted peak: Technical issues

**Solutions:**
```bash
# Re-examine sonication conditions
# Check for nucleosome positioning artifacts
# Validate library preparation protocol
```

### Issue 4: Poor Replicate Consistency

**Diagnosis:**
```bash
# Calculate correlation between replicates
multiBamSummary bins \
    --bamfiles ARF27_rep1.bam ARF27_rep2.bam \
    --outFileName replicate_correlation.npz

plotCorrelation \
    --corData replicate_correlation.npz \
    --plotFile replicate_correlation.png \
    --corMethod pearson
```

## Advanced Quality Metrics

### 1. Signal-to-Noise Ratio Analysis

```bash
# Calculate signal-to-noise using peak summits
computeMatrix reference-point \
    --referencePoint center \
    --beforeRegionStartLength 1000 \
    --afterRegionStartLength 1000 \
    --scoreFileName 05_processed_bams/normalized_bigwigs/ARF27_ChIP_merged_RPM.bw \
    --regionsFileName 07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed \
    --out 08_peak_qc/signal_noise_matrix.gz

# Create signal profile plot
plotProfile \
    --matrixFile 08_peak_qc/signal_noise_matrix.gz \
    --plotFile 08_peak_qc/signal_noise_profile.png \
    --plotTitle "ARF27 Signal Profile at Peak Centers"
```

## Summary Statistics Generation

```bash
# Generate final summary statistics file
echo "Generating final summary statistics..."

cat > 08_peak_qc/final_QC_summary.txt << EOF
ARF27 ChIP-seq Final Quality Summary
====================================
Analysis completed: $(date)

EXPERIMENT OVERVIEW
===================
Samples analyzed: 4 (2 ChIP replicates, 2 Input controls)
Genome: Maize B73-v4 
Peak caller: MACS2
Post-processing: Blacklist filtered

PEAK STATISTICS
===============
Individual replicates:
- Rep1: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak) peaks
- Rep2: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak) peaks

Merged dataset:
- Total peaks: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
- High-confidence peaks (score>100): $(awk '$5 > 100' 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak | wc -l)

QUALITY METRICS
===============
FRiP Scores:
$(cat 08_peak_qc/frip_analysis/frip_scores.txt)

Replicate Reproducibility:
- Overlap percentage: $overlap_rate_rep1%

FILES GENERATED
===============
1. Quality control plots and analyses: 08_peak_qc/
2. Filtered peak files: 07_blacklists/filtered_peaks/
3. Comprehensive QC report: 08_peak_qc/ARF27_ChIP_seq_QC_report.txt
4. Quality checklist: 08_peak_qc/quality_checklist.txt

NEXT STEPS
==========
✓ Quality control completed
→ Ready for Step 9: Differential Peak Analysis

EOF

echo "Final QC summary: 08_peak_qc/final_QC_summary.txt"
```

## Analysis Log Update

```bash
# Update main analysis log
cat >> logs/step8_peak_qc_log.txt << EOF

Step 8: Peak Data Quality Control Results
Date: $(date)

Quality Control Summary:
- FRiP scores calculated for all samples
- Fragment length analysis completed
- Peak enrichment metrics analyzed
- Genomic distribution assessed
- Replicate reproducibility evaluated

Key Findings:
- Total filtered peaks: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
- Best FRiP score: $frip_chip% (merged ChIP)
- Replicate overlap: $overlap_rate_rep1%
- Peak width median: [From analysis] bp

Quality Assessment:
$(cat 08_peak_qc/quality_checklist.txt | grep "OVERALL ASSESSMENT" -A 2)

Files Generated:
- QC plots and analyses: 15+ files
- Comprehensive QC report: 1 file
- Quality checklist: 1 file
- Summary statistics: 2 files

Issues Encountered: [None/List problems and solutions]

Ready for Step 9: Differential Peak Analysis

EOF
```

## Key Takeaways

- Comprehensive QC validates experimental success and data quality
- FRiP scores are the most important single metric for ChIP-seq quality
- Replicate consistency indicates technical reproducibility
- Peak characteristics should match expectations for transcription factors
- Multiple complementary metrics provide robust quality assessment
- Visual inspection of tracks complements statistical metrics

## Next Steps Checklist

Before proceeding to Step 9 (Differential Peak Analysis):

- [ ] FRiP scores calculated and interpreted
- [ ] Fragment length distribution analyzed
- [ ] Peak signal enrichment assessed  
- [ ] Genomic distribution of peaks examined
- [ ] Replicate reproducibility evaluated
- [ ] Peak characteristics (width, shape) analyzed
- [ ] Comprehensive QC report generated
- [ ] Quality checklist completed
- [ ] Any quality issues documented and addressed

Ready for **Step 9: Differential Peak Analysis**? We'll identify condition-specific ARF27 binding sites and quantify changes in binding between experimental conditions using statistical methods designed for ChIP-seq data.
