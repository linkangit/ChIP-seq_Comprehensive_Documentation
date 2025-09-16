# Step 7: Peak Blacklist Filtering for ChIP-seq Data

## Objective
Remove peaks that fall in problematic genomic regions (blacklisted regions) to improve the specificity and reliability of ARF27 binding site predictions by filtering out artifacts and false positives.

## Why Blacklist Filtering is Essential

### What are Blacklisted Regions?
Blacklisted regions are genomic areas known to produce ChIP-seq artifacts due to:
- **Repetitive sequences**: Transposons, tandem repeats, satellite DNA
- **High copy number regions**: Ribosomal RNA genes, histone gene clusters
- **Assembly gaps**: Unresolved regions in reference genome
- **Artifactual regions**: Known to cause false positive signals

### Impact on ChIP-seq Analysis:
- **False positive peaks**: Up to 20-30% of peaks can be artifacts
- **Biased enrichment**: Some regions attract reads non-specifically
- **Downstream analysis errors**: Affect gene annotation and motif discovery
- **Publication concerns**: Reviewers expect blacklist filtering

### Maize-Specific Challenges:
- **High repetitive content**: ~85% of maize genome is repetitive
- **Complex genome structure**: Extensive segmental duplications
- **Centromeric regions**: Large heterochromatic blocks
- **rRNA gene clusters**: Multiple copies causing mapping artifacts

## Understanding Blacklist Sources

### 1. ENCODE Blacklists (Human/Mouse):
- Well-established for model organisms
- Not directly applicable to maize
- Conceptual framework for plant genomes

### 2. Plant-Specific Blacklists:
- Limited availability for crop species
- Need custom generation for maize B73-v4
- Based on repetitive element annotations

### 3. Custom Blacklist Generation:
- High-depth input control analysis
- Repetitive element masking
- Known problematic loci

## Step-by-Step Instructions

### 1. Create Blacklist Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create blacklist directory
mkdir -p 07_blacklists
mkdir -p 07_blacklists/raw_blacklists
mkdir -p 07_blacklists/filtered_peaks
```

### 2. Obtain Repetitive Element Annotations

First, download repetitive element annotations for maize B73-v4:

```bash
# Download repeat annotations from MaizeGDB or RepeatMasker
cd 07_blacklists/raw_blacklists

# Option 1: RepeatMasker output (if available)
# wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.repeatmasker_4.0.9_2019-05-13.gff3

# Option 2: Use TEs from genome annotation
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3

# Extract repetitive elements from GFF3
grep -E "(transposable_element|repetitive_region|tandem_repeat)" Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 > maize_repetitive_elements.gff3

# Convert to BED format
awk '$3=="transposable_element" || $3=="repetitive_region" {print $1"\t"$4-1"\t"$5}' maize_repetitive_elements.gff3 | \
sort -k1,1 -k2,2n > maize_repetitive_elements.bed
```

### 3. Create Custom Blacklist from Input Controls

Generate a blacklist based on high-signal regions in input controls:

```bash
# Create coverage tracks from merged input controls
echo "Creating input control coverage track..."
bamCoverage \
    --bam ../06_peaks/merged_replicates/Input_control_merged.bam \
    --outFileName input_control_coverage.bw \
    --outFileFormat bigwig \
    --binSize 100 \
    --normalizeUsing CPM \
    --effectiveGenomeSize 2100000000 \
    --numberOfProcessors 8

# Identify high-coverage regions in input (potential artifacts)
# Convert BigWig to bedGraph for processing
bigWigToBedGraph input_control_coverage.bw input_control_coverage.bedgraph

# Find regions with >95th percentile coverage
R --vanilla << 'EOF'
# Read coverage data
coverage <- read.table("input_control_coverage.bedgraph", header=FALSE)
colnames(coverage) <- c("chr", "start", "end", "score")

# Calculate 95th percentile threshold
threshold_95 <- quantile(coverage$score, 0.95, na.rm=TRUE)
cat("95th percentile coverage threshold:", threshold_95, "\n")

# Identify high-coverage regions
high_coverage <- coverage[coverage$score > threshold_95, ]

# Merge nearby regions (within 500bp)
if(nrow(high_coverage) > 0) {
    write.table(high_coverage[,1:3], "input_high_coverage_regions.bed", 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
}
EOF

# Merge nearby regions using bedtools
bedtools merge -d 500 -i input_high_coverage_regions.bed > input_blacklist_raw.bed
```

### 4. Create Comprehensive Blacklist

Combine different sources to create a comprehensive blacklist:

```bash
# Combine repetitive elements and input-derived blacklist
echo "Creating comprehensive blacklist..."

# Start with repetitive elements
cp maize_repetitive_elements.bed comprehensive_blacklist_temp.bed

# Add input-derived blacklist if it exists
if [ -f "input_blacklist_raw.bed" ]; then
    cat input_blacklist_raw.bed >> comprehensive_blacklist_temp.bed
fi

# Add known problematic regions (centromeres, rRNA clusters)
# Create manual blacklist for known problematic loci
cat > manual_blacklist.bed << 'EOF'
chr1	147000000	153000000	centromere1
chr2	95000000	100000000	centromere2
chr3	103000000	108000000	centromere3
chr4	106000000	111000000	centromere4
chr5	102000000	107000000	centromere5
chr6	58000000	63000000	centromere6
chr7	58000000	63000000	centromere7
chr8	52000000	57000000	centromere8
chr9	66000000	71000000	centromere9
chr10	44000000	49000000	centromere10
EOF

# Add manual blacklist
cat manual_blacklist.bed >> comprehensive_blacklist_temp.bed

# Sort and merge overlapping regions
sort -k1,1 -k2,2n comprehensive_blacklist_temp.bed | \
bedtools merge -i - > maize_B73_v4_blacklist.bed

# Clean up temporary files
rm comprehensive_blacklist_temp.bed

# Check final blacklist
echo "Final blacklist statistics:"
wc -l maize_B73_v4_blacklist.bed
awk '{sum += $3-$2} END {print "Total blacklisted bases:", sum}' maize_B73_v4_blacklist.bed
```

### 5. Alternative: Download Pre-existing Plant Blacklists

If available, download established plant blacklists:

```bash
# Example: Download from plant ChIP-seq databases (if available)
# wget https://plantregmap.gao-lab.org/download/blacklists/maize_blacklist.bed

# Or use Arabidopsis blacklist as template (for comparative analysis)
# wget https://github.com/Boyle-Lab/Blacklist/raw/master/lists/Athaliana/Athaliana-blacklist.bed
```

### 6. Filter Peaks Using Blacklist

Apply the blacklist to remove peaks in problematic regions:

```bash
# Filter individual replicate peaks
echo "Filtering peaks with blacklist..."

# ARF27 ChIP replicate 1
bedtools intersect \
    -a ../06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak \
    -b raw_blacklists/maize_B73_v4_blacklist.bed \
    -v > ../07_blacklists/filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak

# ARF27 ChIP replicate 2
bedtools intersect \
    -a ../06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak \
    -b raw_blacklists/maize_B73_v4_blacklist.bed \
    -v > ../07_blacklists/filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak

# Merged replicates
bedtools intersect \
    -a ../06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak \
    -b raw_blacklists/maize_B73_v4_blacklist.bed \
    -v > ../07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak

# Filter summits as well
bedtools intersect \
    -a ../06_peaks/merged_replicates/ARF27_ChIP_merged_summits.bed \
    -b raw_blacklists/maize_B73_v4_blacklist.bed \
    -v > ../07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed
```

**Command Explanation:**
- `bedtools intersect`: Find overlapping regions
- `-a`: Query file (peaks)
- `-b`: Target file (blacklist)
- `-v`: Report entries in A that have no overlap with B (filter out)

### 7. Analyze Filtering Impact

Assess how many peaks were removed and their characteristics:

```bash
# Compare peak counts before and after filtering
echo "Blacklist filtering impact analysis:"
echo "Sample\tOriginal_Peaks\tFiltered_Peaks\tRemoved_Peaks\tPercent_Removed"

# Individual replicates
for rep in rep1 rep2
do
    original=$(wc -l < ../06_peaks/individual_replicates/ARF27_ChIP_${rep}_peaks.narrowPeak)
    filtered=$(wc -l < filtered_peaks/ARF27_ChIP_${rep}_peaks_filtered.narrowPeak)
    removed=$((original - filtered))
    percent_removed=$(echo "scale=1; $removed * 100 / $original" | bc)
    echo -e "ARF27_ChIP_${rep}\t${original}\t${filtered}\t${removed}\t${percent_removed}%"
done

# Merged replicates
original=$(wc -l < ../06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak)
filtered=$(wc -l < filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
removed=$((original - filtered))
percent_removed=$(echo "scale=1; $removed * 100 / $original" | bc)
echo -e "ARF27_ChIP_merged\t${original}\t${filtered}\t${removed}\t${percent_removed}%"
```

### 8. Analyze Removed Peaks

Examine characteristics of filtered peaks:

```bash
# Find peaks that were removed by blacklist filtering
bedtools intersect \
    -a ../06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak \
    -b raw_blacklists/maize_B73_v4_blacklist.bed \
    -wa > removed_peaks.narrowPeak

echo "Analysis of removed peaks:"
echo "Number of removed peaks: $(wc -l < removed_peaks.narrowPeak)"

# Analyze removed peak characteristics
R --vanilla << 'EOF'
# Read removed peaks
if(file.exists("removed_peaks.narrowPeak") && file.size("removed_peaks.narrowPeak") > 0) {
    removed <- read.table("removed_peaks.narrowPeak", sep="\t")
    colnames(removed) <- c("chr", "start", "end", "name", "score", "strand", 
                          "signalValue", "pValue", "qValue", "peak")
    
    # Read filtered (kept) peaks
    filtered <- read.table("filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak", sep="\t")
    colnames(filtered) <- c("chr", "start", "end", "name", "score", "strand", 
                           "signalValue", "pValue", "qValue", "peak")
    
    cat("Removed peaks characteristics:\n")
    cat("Mean score (removed):", round(mean(removed$score), 1), "\n")
    cat("Mean score (kept):", round(mean(filtered$score), 1), "\n")
    cat("Median length (removed):", median(removed$end - removed$start), "bp\n")
    cat("Median length (kept):", median(filtered$end - filtered$start), "bp\n")
    
    # Score comparison
    png("../07_blacklists/removed_vs_kept_scores.png", width=800, height=600)
    boxplot(list(Removed=removed$score, Kept=filtered$score),
            main="Peak Score Comparison: Removed vs Kept",
            ylab="MACS2 Score", col=c("lightcoral", "lightblue"))
    dev.off()
    
    cat("Score comparison plot saved\n")
} else {
    cat("No peaks were removed by blacklist filtering\n")
}
EOF
```

## Quality Control and Validation

### 1. Validate Blacklist Coverage

```bash
# Calculate genome coverage of blacklist
total_blacklist_bp=$(awk '{sum += $3-$2} END {print sum}' raw_blacklists/maize_B73_v4_blacklist.bed)
genome_size=2100000000
blacklist_percent=$(echo "scale=2; $total_blacklist_bp * 100 / $genome_size" | bc)

echo "Blacklist coverage analysis:"
echo "Total blacklisted bases: $total_blacklist_bp"
echo "Genome size: $genome_size"
echo "Blacklist coverage: $blacklist_percent%"

# Typical values: 5-15% for plant genomes with high repeat content
```

### 2. Chromosome-wise Analysis

```bash
# Analyze blacklist distribution across chromosomes
echo "Blacklist distribution by chromosome:"
echo "Chromosome\tBlacklisted_Regions\tBlacklisted_Bases\tPercent_of_Chr"

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
do
    regions=$(awk -v chr="$chr" '$1==chr' raw_blacklists/maize_B73_v4_blacklist.bed | wc -l)
    bases=$(awk -v chr="$chr" '$1==chr {sum += $3-$2} END {print sum+0}' raw_blacklists/maize_B73_v4_blacklist.bed)
    
    # Approximate chromosome lengths for maize B73-v4
    case $chr in
        chr1) chr_length=307041717 ;;
        chr2) chr_length=244442276 ;;
        chr3) chr_length=235667834 ;;
        chr4) chr_length=246994605 ;;
        chr5) chr_length=223902240 ;;
        chr6) chr_length=174033170 ;;
        chr7) chr_length=182381542 ;;
        chr8) chr_length=175161014 ;;
        chr9) chr_length=159769782 ;;
        chr10) chr_length=150982314 ;;
    esac
    
    percent=$(echo "scale=2; $bases * 100 / $chr_length" | bc)
    echo -e "${chr}\t${regions}\t${bases}\t${percent}%"
done
```

### 3. Peak Quality Assessment Post-Filtering

```bash
# Assess if filtering improved peak quality
echo "Peak quality assessment after blacklist filtering:"

# Compare replicate overlap before and after filtering
echo "Replicate overlap analysis:"

# Before filtering
bedtools intersect \
    -a ../06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak \
    -b ../06_peaks/individual_replicates/ARF27_ChIP_rep2_peaks.narrowPeak \
    -wa | wc -l > temp_overlap_before.txt

rep1_before=$(wc -l < ../06_peaks/individual_replicates/ARF27_ChIP_rep1_peaks.narrowPeak)
overlap_before=$(cat temp_overlap_before.txt)
percent_before=$(echo "scale=1; $overlap_before * 100 / $rep1_before" | bc)

# After filtering
bedtools intersect \
    -a filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak \
    -b filtered_peaks/ARF27_ChIP_rep2_peaks_filtered.narrowPeak \
    -wa | wc -l > temp_overlap_after.txt

rep1_after=$(wc -l < filtered_peaks/ARF27_ChIP_rep1_peaks_filtered.narrowPeak)
overlap_after=$(cat temp_overlap_after.txt)
percent_after=$(echo "scale=1; $overlap_after * 100 / $rep1_after" | bc)

echo "Before filtering: $overlap_before/$rep1_before ($percent_before%) overlap"
echo "After filtering: $overlap_after/$rep1_after ($percent_after%) overlap"

# Clean up
rm temp_overlap_before.txt temp_overlap_after.txt
```

## Advanced Blacklist Customization

### 1. Create Signal-to-Noise Based Blacklist

For more sophisticated filtering:

```bash
# Calculate signal-to-noise ratio for each peak
R --vanilla << 'EOF'
# This requires the MACS2 fold-enrichment track
# Read peak data
peaks <- read.table("../06_peaks/merged_replicates/ARF27_ChIP_merged_peaks.narrowPeak", sep="\t")
colnames(peaks) <- c("chr", "start", "end", "name", "score", "strand", 
                    "signalValue", "pValue", "qValue", "peak")

# Identify peaks with low signal-to-noise (fold-change < 2)
low_quality <- peaks[peaks$signalValue < 2, ]

if(nrow(low_quality) > 0) {
    write.table(low_quality[,1:3], "low_quality_peaks.bed", 
                sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    cat("Identified", nrow(low_quality), "low quality peaks\n")
}
EOF

# Add low-quality peaks to blacklist (optional)
if [ -f "low_quality_peaks.bed" ]; then
    cat raw_blacklists/maize_B73_v4_blacklist.bed low_quality_peaks.bed | \
    sort -k1,1 -k2,2n | \
    bedtools merge -i - > raw_blacklists/maize_B73_v4_blacklist_extended.bed
fi
```

### 2. Gene-based Blacklist Filtering

Exclude peaks in specific gene families (e.g., rRNA, tRNA):

```bash
# Download gene annotation if not already available
# Extract rRNA and tRNA genes
awk '$3=="rRNA" || $3=="tRNA"' ../reference_genome/maize_annotation.gff3 | \
awk '{print $1"\t"$4-1"\t"$5}' > rRNA_tRNA_genes.bed

# Optionally add to blacklist
# cat raw_blacklists/maize_B73_v4_blacklist.bed rRNA_tRNA_genes.bed | \
# sort -k1,1 -k2,2n | bedtools merge -i - > enhanced_blacklist.bed
```

## Common Issues and Troubleshooting

### Issue 1: Too Many Peaks Removed (>50%)

**Possible Causes:**
- Overly aggressive blacklist
- Poor ChIP specificity
- Wrong blacklist file

**Solutions:**
```bash
# Check blacklist size
awk '{sum += $3-$2} END {print "Blacklist covers", sum/1000000, "Mb"}' raw_blacklists/maize_B73_v4_blacklist.bed

# Use more conservative blacklist (only centromeres and known artifacts)
grep -E "(centromere|artifact)" raw_blacklists/maize_B73_v4_blacklist.bed > conservative_blacklist.bed

# Re-filter with conservative blacklist
bedtools intersect -a peaks.narrowPeak -b conservative_blacklist.bed -v > peaks_conservative_filtered.narrowPeak
```

### Issue 2: Very Few Peaks Removed (<5%)

**Possible Causes:**
- Insufficient blacklist
- High-quality ChIP with good specificity
- Missing repetitive element annotations

**Assessment:**
```bash
# Check if blacklist is comprehensive
echo "Blacklist validation:"
echo "Total regions: $(wc -l < raw_blacklists/maize_B73_v4_blacklist.bed)"
echo "Coverage: $blacklist_percent% of genome"

# If coverage is <5%, blacklist may be insufficient
```

### Issue 3: Blacklist File Format Issues

**Common problems:**
- Wrong coordinate system (0-based vs 1-based)
- Unsorted file
- Extra columns

**Solutions:**
```bash
# Ensure proper BED format (0-based, sorted)
awk 'NF>=3 {print $1"\t"$2"\t"$3}' input_blacklist.bed | \
sort -k1,1 -k2,2n > fixed_blacklist.bed

# Check for coordinate issues
awk '$2 >= $3 {print "Bad coordinates:", $0}' fixed_blacklist.bed
```

### Issue 4: Memory Issues with Large Blacklists

**Solutions:**
```bash
# Process chromosomes individually
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
do
    # Filter peaks for each chromosome
    awk -v chr="$chr" '$1==chr' peaks.narrowPeak | \
    bedtools intersect -a - -b <(awk -v chr="$chr" '$1==chr' blacklist.bed) -v \
    >> peaks_filtered.narrowPeak
done
```

## Validation and Best Practices

### 1. Cross-validate with Known Biology

```bash
# Check if known ARF binding sites are retained
# Example: If you have validated ARF27 targets
# bedtools intersect -a known_ARF27_sites.bed -b filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak -wa
```

### 2. Compare with Literature

```bash
# Extract peaks near auxin-responsive genes
# awk '$3=="gene" && $9 ~ /auxin/' annotation.gff3 | bedtools window -a filtered_peaks.narrowPeak -b - -w 5000
```

### 3. Document Filtering Decisions

```bash
# Create filtering report
cat > ../07_blacklists/filtering_report.txt << EOF
Blacklist Filtering Report
Date: $(date)

Blacklist Sources:
- Repetitive elements: $(wc -l < raw_blacklists/maize_repetitive_elements.bed) regions
- Input control artifacts: $(wc -l < raw_blacklists/input_blacklist_raw.bed 2>/dev/null || echo "0") regions
- Manual annotations: $(wc -l < raw_blacklists/manual_blacklist.bed) regions

Final Blacklist:
- Total regions: $(wc -l < raw_blacklists/maize_B73_v4_blacklist.bed)
- Genome coverage: $blacklist_percent%

Filtering Impact:
- Original peaks: $original
- Filtered peaks: $filtered
- Removed peaks: $removed ($percent_removed%)

Quality Assessment:
- Replicate overlap improved: [Yes/No]
- Peak scores maintained: [Yes/No]
- Biological relevance preserved: [Yes/No]

Files Generated:
- Blacklist: raw_blacklists/maize_B73_v4_blacklist.bed
- Filtered peaks: filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak
- Filtered summits: filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed
EOF
```

## Analysis Log Update

```bash
# Update analysis log
cat >> logs/step7_blacklist_filtering_log.txt << EOF

Step 7: Blacklist Filtering Results
Date: $(date)

Blacklist Creation:
- Sources used: Repetitive elements, input artifacts, manual curation
- Total blacklisted regions: $(wc -l < 07_blacklists/raw_blacklists/maize_B73_v4_blacklist.bed)
- Genome coverage: $blacklist_percent%

Filtering Results:
$(echo "Sample\tOriginal\tFiltered\tRemoved\tPercent_Removed")
$(grep -v "Sample" filtering_impact.txt)

Quality Impact:
- Replicate overlap before: $percent_before%
- Replicate overlap after: $percent_after%
- Peak quality: [Improved/Maintained/Degraded]

Files Generated:
- Comprehensive blacklist: 1 file
- Filtered peak files: 3 files
- Quality control plots: Multiple files

Issues Encountered: [None/List problems and solutions]

Ready for Step 8: Peak Data Quality Control

EOF
```

## Key Takeaways

- Blacklist filtering is essential for removing false positive peaks
- Custom blacklists are often needed for non-model organisms like maize
- Typical removal rate: 10-30% of peaks in repetitive genomes
- Quality should improve (better replicate overlap) after filtering
- Document all filtering decisions for reproducibility

## Next Steps Checklist

Before proceeding to Step 8 (Peak Data Quality Control):

- [ ] Comprehensive blacklist created and validated
- [ ] All peak files successfully filtered
- [ ] Filtering impact assessed and documented
- [ ] Peak quality metrics compared before/after filtering
- [ ] Removed peaks analyzed for characteristics
- [ ] Final filtered peak files organized and verified

Ready for **Step 8: Peak Data Quality Control**? We'll perform comprehensive quality assessment of the filtered peaks including fragment length analysis, signal enrichment metrics, and peak annotation to ensure our ARF27 binding sites are of high quality and biologically relevant.
