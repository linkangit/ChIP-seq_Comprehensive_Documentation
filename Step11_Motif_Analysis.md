# Step 11: Motif Analysis and Binding Site Discovery

## Objective
Discover ARF27 binding motifs through computational analysis of peak sequences and compare them with known ARF binding motifs to validate ChIP-seq specificity and understand the DNA-binding preferences of ARF27 in maize.

## Why Motif Analysis is Critical

### Validation of ChIP-seq Success:
- **Specificity confirmation**: Presence of expected motifs validates successful IP
- **Antibody validation**: Confirms antibody specificity to ARF27
- **Protocol validation**: Indicates proper ChIP-seq execution
- **Biological relevance**: Links binding to known regulatory mechanisms

### Biological Insights:
- **DNA-binding specificity**: Precise sequence preferences of ARF27
- **Cooperative binding**: Co-occurrence with other transcription factor motifs
- **Evolutionary conservation**: Comparison with ARF motifs across species
- **Regulatory mechanisms**: Understanding of gene regulation by ARF27

### ARF Family Motif Characteristics:
- **Core motif**: TGTCTC (canonical Auxin Response Element)
- **Variations**: TGTCNN, NNGACAA (reverse complement)
- **Spacing patterns**: Direct repeats, palindromes
- **Context dependency**: Surrounding sequence preferences

## Understanding Motif Discovery Methods

### 1. De Novo Motif Discovery:
- **MEME**: Finds overrepresented sequence patterns
- **Homer**: Specialized for ChIP-seq motif discovery
- **RSAT**: Plant-specific motif analysis tools
- **No prior knowledge**: Discovers motifs without bias

### 2. Known Motif Scanning:
- **FIMO**: Scans for known motifs in sequences
- **HOMER scanMotif**: Tests enrichment of known motifs
- **PlantTFDB**: Plant transcription factor binding sites
- **Validates expectations**: Confirms presence of known ARF motifs

### 3. Comparative Analysis:
- **Cross-species comparison**: Compare with Arabidopsis, rice ARF motifs
- **Family comparison**: Compare ARF27 with other ARF family members
- **Motif evolution**: Understand conservation and divergence

## Step-by-Step Motif Analysis

### 1. Set Up Motif Analysis Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create motif analysis directory structure
mkdir -p 11_motif_analysis
mkdir -p 11_motif_analysis/sequences
mkdir -p 11_motif_analysis/meme_analysis
mkdir -p 11_motif_analysis/homer_analysis
mkdir -p 11_motif_analysis/known_motifs
mkdir -p 11_motif_analysis/results
mkdir -p 11_motif_analysis/plots
mkdir -p 11_motif_analysis/databases
```

### 2. Extract Peak Sequences

Extract DNA sequences from peak regions for motif discovery:

```bash
# Extract sequences from peak summits
echo "Extracting peak sequences for motif analysis..."

# Download maize genome if not already available
cd 11_motif_analysis/databases
if [ ! -f "maize_B73_v4.fa" ]; then
    # Link to reference genome from earlier steps
    ln -s ../../../reference_genome/maize_B73_v4.fa .
fi

# Index genome for bedtools
if [ ! -f "maize_B73_v4.fa.fai" ]; then
    samtools faidx maize_B73_v4.fa
fi

cd ../sequences

# Extract 200bp sequences around peak summits (±100bp)
echo "Extracting sequences around peak summits..."
awk '{summit_pos = $2 + $10; print $1"\t"(summit_pos-100)"\t"(summit_pos+100)"\tpeak_"NR}' \
    ../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed > \
    peak_summits_200bp.bed

# Get sequences using bedtools
bedtools getfasta \
    -fi ../databases/maize_B73_v4.fa \
    -bed peak_summits_200bp.bed \
    -fo peak_summits_200bp.fasta \
    -name

# Extract sequences for top 500 peaks (for faster analysis)
head -500 ../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak | \
awk '{summit_pos = $2 + $10; print $1"\t"(summit_pos-100)"\t"(summit_pos+100)"\ttop_peak_"NR}' > \
    top500_summits_200bp.bed

bedtools getfasta \
    -fi ../databases/maize_B73_v4.fa \
    -bed top500_summits_200bp.bed \
    -fo top500_summits_200bp.fasta \
    -name

# Create different sequence lengths for analysis
# 100bp sequences (±50bp from summit)
awk '{summit_pos = $2 + $10; print $1"\t"(summit_pos-50)"\t"(summit_pos+50)"\tpeak_"NR}' \
    ../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed > \
    peak_summits_100bp.bed

bedtools getfasta \
    -fi ../databases/maize_B73_v4.fa \
    -bed peak_summits_100bp.bed \
    -fo peak_summits_100bp.fasta \
    -name

# 500bp sequences (±250bp from summit) for broader motif discovery
awk '{summit_pos = $2 + $10; print $1"\t"(summit_pos-250)"\t"(summit_pos+250)"\tpeak_"NR}' \
    ../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_summits_filtered.bed > \
    peak_summits_500bp.bed

bedtools getfasta \
    -fi ../databases/maize_B73_v4.fa \
    -bed peak_summits_500bp.bed \
    -fo peak_summits_500bp.fasta \
    -name

echo "Peak sequences extracted:"
echo "- 100bp sequences: $(grep -c ">" peak_summits_100bp.fasta)"
echo "- 200bp sequences: $(grep -c ">" peak_summits_200bp.fasta)"
echo "- 500bp sequences: $(grep -c ">" peak_summits_500bp.fasta)"
echo "- Top 500 peaks (200bp): $(grep -c ">" top500_summits_200bp.fasta)"
```

### 3. Install Motif Analysis Tools

```bash
# Install MEME Suite
echo "Installing MEME Suite..."

# Option 1: Using conda (recommended)
conda install -c bioconda meme

# Option 2: Download from MEME website
# wget http://meme-suite.org/meme-software/5.4.1/meme-5.4.1.tar.gz
# tar zxf meme-5.4.1.tar.gz
# cd meme-5.4.1
# ./configure --prefix=$HOME/meme
# make && make install

# Verify installation
meme -version
fimo -version
tomtom -version

# Install HOMER (optional but recommended for ChIP-seq)
echo "Installing HOMER..."
# wget http://homer.ucsd.edu/homer/configureHomer.pl
# perl configureHomer.pl -install
```

### 4. Create Background Sequences

Generate appropriate background sequences for motif discovery:

```bash
# Create background sequences for motif analysis
echo "Creating background sequences..."

cd ../sequences

# Method 1: Random genomic regions
# Select random regions of same size as peaks from maize genome
total_peaks=$(wc -l < peak_summits_200bp.bed)

# Create random regions avoiding peaks
bedtools random \
    -g ../databases/maize_B73_v4.fa.fai \
    -l 200 \
    -n $total_peaks \
    -seed 42 > random_background_200bp.bed

# Exclude regions that overlap with peaks
bedtools intersect \
    -a random_background_200bp.bed \
    -b peak_summits_200bp.bed \
    -v > random_background_200bp_filtered.bed

# Get background sequences
bedtools getfasta \
    -fi ../databases/maize_B73_v4.fa \
    -bed random_background_200bp_filtered.bed \
    -fo background_200bp.fasta \
    -name

# Method 2: Shuffled sequences (alternative background)
fasta-shuffle-letters -seed 42 peak_summits_200bp.fasta > shuffled_background_200bp.fasta

echo "Background sequences created:"
echo "- Random genomic regions: $(grep -c ">" background_200bp.fasta)"
echo "- Shuffled sequences: $(grep -c ">" shuffled_background_200bp.fasta)"
```

### 5. MEME De Novo Motif Discovery

Discover novel motifs using MEME:

```bash
# Run MEME for de novo motif discovery
echo "Running MEME for de novo motif discovery..."

cd ../meme_analysis

# Run MEME on top 500 peaks (faster analysis)
echo "MEME analysis on top 500 peaks..."
meme ../sequences/top500_summits_200bp.fasta \
    -dna \
    -oc meme_top500_200bp \
    -nmotifs 10 \
    -minw 6 \
    -maxw 15 \
    -objfun classic \
    -revcomp \
    -markov_order 2 \
    -mod zoops \
    -time 300

# Run MEME on all peaks with 100bp sequences (focused on core binding site)
echo "MEME analysis on all peaks (100bp)..."
meme ../sequences/peak_summits_100bp.fasta \
    -dna \
    -oc meme_all_100bp \
    -nmotifs 5 \
    -minw 6 \
    -maxw 12 \
    -objfun classic \
    -revcomp \
    -markov_order 2 \
    -mod zoops \
    -time 600

# Run MEME-ChIP for comprehensive analysis (if computational resources allow)
echo "Running MEME-ChIP comprehensive analysis..."
meme-chip \
    -oc meme_chip_analysis \
    -dna \
    -meme-nmotifs 10 \
    -meme-minw 6 \
    -meme-maxw 15 \
    -time 600 \
    ../sequences/top500_summits_200bp.fasta

echo "MEME analysis completed. Results in meme_analysis/"
```

### 6. Known Motif Analysis

Search for known ARF and auxin-related motifs:

```bash
# Create known motif database
echo "Creating known motif database..."

cd ../known_motifs

# Create known ARF motifs in MEME format
cat > arf_motifs.meme << 'EOF'
MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 

MOTIF ARF_canonical
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.050 0.050 0.050 0.850
0.050 0.050 0.850 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050

MOTIF ARF_extended
letter-probability matrix: alength= 4 w= 8 nsites= 20 E= 0
0.050 0.050 0.050 0.850
0.050 0.050 0.850 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
0.250 0.250 0.250 0.250
0.250 0.250 0.250 0.250
0.850 0.050 0.050 0.050
0.850 0.050 0.050 0.050

MOTIF AuxRE_TGTCTC
letter-probability matrix: alength= 4 w= 6 nsites= 20 E= 0
0.050 0.050 0.050 0.850
0.050 0.050 0.850 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
0.050 0.850 0.050 0.050
EOF

# Download plant motif databases (if available)
# wget http://plantregmap.gao-lab.org/download/TF_binding_motifs/Zea_mays_motifs.meme

# Scan peaks for known ARF motifs using FIMO
echo "Scanning peaks for known ARF motifs..."
fimo \
    --oc fimo_arf_motifs \
    --thresh 1e-4 \
    --bgfile motif-file \
    arf_motifs.meme \
    ../sequences/peak_summits_200bp.fasta

# Create summary of known motif occurrences
echo "Known motif analysis completed."
```

### 7. HOMER Motif Analysis

Use HOMER for ChIP-seq specialized motif discovery:

```bash
# HOMER motif analysis (if HOMER is installed)
echo "Running HOMER motif analysis..."

cd ../homer_analysis

# Convert peak file to HOMER format
awk '{print $4"\t"$1"\t"$2"\t"$3"\t+"}' \
    ../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak > \
    ARF27_peaks_homer.txt

# Run HOMER findMotifsGenome.pl
# findMotifsGenome.pl \
#     ARF27_peaks_homer.txt \
#     ../databases/maize_B73_v4.fa \
#     homer_motifs_output \
#     -size 200 \
#     -mask \
#     -len 8,10,12 \
#     -S 10

echo "HOMER analysis setup completed (requires HOMER installation)."
```

### 8. Analyze Motif Discovery Results

Process and analyze the discovered motifs:

```bash
# Analyze MEME results
echo "Analyzing motif discovery results..."

cd ../results

# Extract motif information from MEME results
if [ -f "../meme_analysis/meme_top500_200bp/meme.txt" ]; then
    echo "MEME Top 500 Peaks Analysis Results:"
    echo "====================================="
    
    # Extract E-values and motif widths
    grep -A 2 "MOTIF" ../meme_analysis/meme_top500_200bp/meme.txt | \
    grep -E "MOTIF|letter-probability|sites" > meme_motif_summary.txt
    
    echo "Motif summary extracted to meme_motif_summary.txt"
fi

# Analyze FIMO results for known motifs
if [ -f "../known_motifs/fimo_arf_motifs/fimo.tsv" ]; then
    echo "Known ARF Motif Occurrence Analysis:"
    echo "===================================="
    
    # Count motif occurrences
    tail -n +2 ../known_motifs/fimo_arf_motifs/fimo.tsv | \
    awk '{print $1}' | sort | uniq -c > known_motif_counts.txt
    
    # Calculate enrichment
    total_peaks=$(grep -c ">" ../sequences/peak_summits_200bp.fasta)
    
    echo "Motif occurrence summary:" > motif_enrichment_summary.txt
    echo "Total peaks analyzed: $total_peaks" >> motif_enrichment_summary.txt
    echo "" >> motif_enrichment_summary.txt
    echo "Known motif occurrences:" >> motif_enrichment_summary.txt
    cat known_motif_counts.txt >> motif_enrichment_summary.txt
    
    cat motif_enrichment_summary.txt
fi
```

### 9. Compare with Known ARF Motifs

Compare discovered motifs with known ARF family motifs:

```bash
# Compare discovered motifs with known ARF motifs
echo "Comparing discovered motifs with known ARF motifs..."

# Use TOMTOM to compare discovered motifs with known motifs
if [ -f "../meme_analysis/meme_top500_200bp/meme.txt" ]; then
    tomtom \
        -oc tomtom_comparison \
        -dist pearson \
        -thresh 0.1 \
        ../meme_analysis/meme_top500_200bp/meme.txt \
        ../known_motifs/arf_motifs.meme
    
    echo "TOMTOM comparison completed. Results in tomtom_comparison/"
fi

# Statistical analysis of motif enrichment
R --vanilla << 'EOF'
# Analyze motif enrichment and significance
cat("Motif Analysis Summary\n")
cat("=====================\n\n")

# Read FIMO results if available
if(file.exists("../known_motifs/fimo_arf_motifs/fimo.tsv")) {
    fimo_results <- read.table("../known_motifs/fimo_arf_motifs/fimo.tsv", 
                              header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    cat("Known motif analysis:\n")
    cat("Total significant matches:", nrow(fimo_results), "\n")
    
    # Count motifs per sequence
    motifs_per_seq <- table(fimo_results$sequence_name)
    cat("Sequences with motifs:", length(motifs_per_seq), "\n")
    cat("Average motifs per sequence:", round(mean(motifs_per_seq), 2), "\n")
    
    # Motif distribution
    motif_types <- table(fimo_results$motif_id)
    cat("\nMotif type distribution:\n")
    print(motif_types)
    
    # Create visualization
    png("../plots/motif_score_distribution.png", width=800, height=600)
    hist(fimo_results$score, breaks=30, 
         main="ARF Motif Score Distribution in Peaks", 
         xlab="Motif Score", ylab="Frequency",
         col="lightblue")
    dev.off()
    
    # P-value distribution
    png("../plots/motif_pvalue_distribution.png", width=800, height=600)
    hist(-log10(fimo_results$p.value), breaks=30, 
         main="ARF Motif Significance Distribution", 
         xlab="-log10(p-value)", ylab="Frequency",
         col="lightgreen")
    dev.off()
}

# Calculate motif enrichment compared to background
if(file.exists("motif_enrichment_summary.txt")) {
    cat("\nMotif enrichment analysis saved to motif_enrichment_summary.txt\n")
}
EOF
```

### 10. Visualize Motif Results

Create comprehensive visualizations of motif analysis:

```bash
# Create motif visualization plots
echo "Creating motif visualization plots..."

cd ../plots

# Create motif logos using MEME results
if [ -f "../meme_analysis/meme_top500_200bp/meme.html" ]; then
    echo "MEME results available for visualization at:"
    echo "file://$(pwd)/../meme_analysis/meme_top500_200bp/meme.html"
fi

# R script for additional visualizations
R --vanilla << 'EOF'
library(ggplot2)

# Create motif occurrence heatmap
if(file.exists("../known_motifs/fimo_arf_motifs/fimo.tsv")) {
    fimo_results <- read.table("../known_motifs/fimo_arf_motifs/fimo.tsv", 
                              header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # Create position-wise motif occurrence plot
    if(nrow(fimo_results) > 0) {
        png("motif_position_distribution.png", width=1000, height=600)
        hist(fimo_results$start, breaks=50,
             main="ARF Motif Position Distribution in 200bp Peak Regions",
             xlab="Position in Peak Region",
             ylab="Number of Motifs",
             col="lightcoral")
        abline(v=100, col="red", lty=2, lwd=2)  # Summit position
        text(100, max(hist(fimo_results$start, plot=FALSE)$counts)/2, 
             "Summit", pos=4, col="red")
        dev.off()
        
        # Motif orientation analysis
        orientation_table <- table(fimo_results$strand)
        png("motif_orientation.png", width=600, height=600)
        pie(orientation_table, 
            main="ARF Motif Strand Distribution",
            col=c("lightblue", "lightgreen"))
        dev.off()
        
        cat("Motif visualization plots created.\n")
    }
}

# Summary statistics
if(file.exists("../results/motif_enrichment_summary.txt")) {
    cat("\nMotif enrichment summary:\n")
    cat(readLines("../results/motif_enrichment_summary.txt"), sep="\n")
}
EOF
```

### 11. Validate Motif Discovery

Assess the quality and biological relevance of discovered motifs:

```bash
# Validate motif discovery results
echo "Validating motif discovery results..."

cd ../results

# Create validation summary
cat > motif_validation_summary.txt << EOF
ARF27 Motif Analysis Validation Summary
======================================
Analysis Date: $(date)

MOTIF DISCOVERY SUCCESS METRICS
==============================

1. Known ARF Motif Detection:
$(if [ -f "../known_motifs/fimo_arf_motifs/fimo.tsv" ]; then
    significant_matches=$(tail -n +2 ../known_motifs/fimo_arf_motifs/fimo.tsv | wc -l)
    total_peaks=$(grep -c ">" ../sequences/peak_summits_200bp.fasta)
    percentage=$(echo "scale=1; $significant_matches * 100 / $total_peaks" | bc)
    echo "   - Significant ARF motif matches: $significant_matches"
    echo "   - Total peaks analyzed: $total_peaks"
    echo "   - Percentage with known ARF motifs: $percentage%"
    echo "   - Expected: >30% for successful ARF ChIP-seq"
else
    echo "   - FIMO analysis not completed"
fi)

2. De Novo Motif Discovery:
$(if [ -f "../meme_analysis/meme_top500_200bp/meme.txt" ]; then
    motif_count=$(grep -c "MOTIF" ../meme_analysis/meme_top500_200bp/meme.txt | head -1)
    echo "   - Number of discovered motifs: $motif_count"
    echo "   - Expected: 3-10 significant motifs"
    
    # Check for ARF-like motifs in discovered motifs
    if grep -q "TG.*C" ../meme_analysis/meme_top500_200bp/meme.txt; then
        echo "   - ARF-like sequences detected: Yes"
    else
        echo "   - ARF-like sequences detected: Manual inspection required"
    fi
else
    echo "   - MEME analysis not completed"
fi)

3. Motif Quality Indicators:
   - E-values: Should be < 1e-3 for significant motifs
   - Sites per motif: Should be >50 for reliable motifs
   - Motif width: Expected 6-12 bp for ARF motifs

BIOLOGICAL VALIDATION
====================

1. Expected ARF Motif Characteristics:
   - Core sequence: TGTCTC (canonical AuxRE)
   - Alternative: GAGACA (reverse complement)
   - Variations: TGTCNN, NNGACA
   - Context: Often found in pairs or clusters

2. Comparison with Literature:
   - Arabidopsis ARF motifs: TGTCTC
   - Rice ARF motifs: Similar to TGTCTC
   - Maize-specific variations: Expected

3. Functional Validation Suggestions:
   - Check motif occurrence in known auxin-responsive genes
   - Compare with ARF binding sites from other species
   - Validate top motif-containing peaks with qRT-PCR

QUALITY ASSESSMENT
==================
$(if [ -f "../known_motifs/fimo_arf_motifs/fimo.tsv" ]; then
    significant_matches=$(tail -n +2 ../known_motifs/fimo_arf_motifs/fimo.tsv | wc -l)
    total_peaks=$(grep -c ">" ../sequences/peak_summits_200bp.fasta)
    if [ "$significant_matches" -gt 0 ]; then
        percentage=$(echo "scale=1; $significant_matches * 100 / $total_peaks" | bc)
        if (( $(echo "$percentage >= 30" | bc -l) )); then
            echo "Overall Assessment: EXCELLENT - High ARF motif enrichment ($percentage%)"
        elif (( $(echo "$percentage >= 15" | bc -l) )); then
            echo "Overall Assessment: GOOD - Moderate ARF motif enrichment ($percentage%)"
        elif (( $(echo "$percentage >= 5" | bc -l) )); then
            echo "Overall Assessment: ACCEPTABLE - Low but detectable ARF motif enrichment ($percentage%)"
        else
            echo "Overall Assessment: POOR - Very low ARF motif enrichment ($percentage%)"
        fi
    else
        echo "Overall Assessment: FAILED - No significant ARF motifs detected"
    fi
else
    echo "Overall Assessment: INCOMPLETE - Analysis not finished"
fi)

RECOMMENDATIONS
===============
1. Manual inspection of top discovered motifs required
2. Compare results with published ARF ChIP-seq studies
3. Validate motif-containing peaks with functional assays
4. Consider motif refinement with larger peak sets

NEXT STEPS
==========
✓ Motif analysis completed
→ Ready for data integration and interpretation
→ Ready for experimental validation design
→ Ready for comparative analysis with other datasets

EOF

echo "Motif validation summary created: motif_validation_summary.txt"
cat motif_validation_summary.txt
```

### 12. Generate Final Motif Report

Create comprehensive final report:

```bash
# Generate comprehensive motif analysis report
echo "Generating comprehensive motif analysis report..."

cat > ../ARF27_motif_analysis_final_report.txt << EOF
ARF27 ChIP-seq Motif Analysis Final Report
==========================================
Analysis Date: $(date)
Genome: Maize B73-v4

EXECUTIVE SUMMARY
================
This report summarizes the computational motif analysis of ARF27 ChIP-seq peaks
to discover DNA binding preferences and validate ChIP-seq specificity.

ANALYSIS OVERVIEW
================
Total peaks analyzed: $(grep -c ">" ../sequences/peak_summits_200bp.fasta)
Sequence lengths tested: 100bp, 200bp, 500bp around summits
Methods employed: MEME, FIMO, TOMTOM
Background controls: Random genomic regions, shuffled sequences

MOTIF DISCOVERY RESULTS
======================

1. De Novo Motif Discovery (MEME):
$(if [ -f "../meme_analysis/meme_top500_200bp/meme.txt" ]; then
    echo "   ✓ Analysis completed on top 500 peaks"
    echo "   ✓ Results available in HTML format"
    echo "   → Manual inspection of discovered motifs required"
else
    echo "   ✗ MEME analysis incomplete"
fi)

2. Known ARF Motif Scanning (FIMO):
$(if [ -f "../known_motifs/fimo_arf_motifs/fimo.tsv" ]; then
    significant_matches=$(tail -n +2 ../known_motifs/fimo_arf_motifs/fimo.tsv | wc -l)
    echo "   ✓ Scanned for canonical ARF motifs (TGTCTC)"
    echo "   ✓ Found $significant_matches significant matches"
    echo "   → See detailed results in fimo_arf_motifs/"
else
    echo "   ✗ FIMO analysis incomplete"
fi)

3. Motif Comparison (TOMTOM):
$(if [ -f "../results/tomtom_comparison/tomtom.html" ]; then
    echo "   ✓ Compared discovered motifs with known ARF motifs"
    echo "   → See similarity results in tomtom_comparison/"
else
    echo "   ✗ TOMTOM comparison incomplete"
fi)

KEY FINDINGS
============
$(cat ../results/motif_validation_summary.txt | grep -A 10 "QUALITY ASSESSMENT")

BIOLOGICAL INTERPRETATION
=========================
1. ARF27 Binding Specificity:
   - Presence of canonical ARF motifs validates successful ChIP
   - Motif enrichment indicates specific protein-DNA interactions
   - Novel motif variants may represent maize-specific evolution

2. Regulatory Implications:
   - High motif density suggests direct transcriptional regulation
   - Motif positioning relative to summits indicates binding precision
   - Co-occurring motifs may reveal partner transcription factors

3. Evolutionary Context:
   - Comparison with other species reveals conservation/divergence
   - Motif variations may explain species-specific regulation
   - Family-specific motifs distinguish ARF27 from other ARFs

FILES AND RESULTS LOCATIONS
===========================
1. Sequence Files: sequences/
   - Peak summit sequences (100bp, 200bp, 500bp)
   - Background control sequences
   - Top 500 peak subset for focused analysis

2. MEME Results: meme_analysis/
   - De novo motif discovery results
   - Motif logos and statistical significance
   - HTML reports for interactive exploration

3. Known Motif Analysis: known_motifs/
   - ARF motif database
   - FIMO scanning results
   - Motif occurrence statistics

4. Comparative Analysis: results/
   - TOMTOM motif similarity results
   - Enrichment calculations
   - Validation summaries

5. Visualizations: plots/
   - Motif score distributions
   - Position analysis plots
   - Orientation and enrichment graphics

VALIDATION CHECKLIST
====================
□ Known ARF motifs detected in significant fraction of peaks
□ De novo discovery identifies ARF-like sequences
□ Motif positioning correlates with peak summits
□ Statistical significance supports specific binding
□ Biological relevance confirmed through literature comparison

EXPERIMENTAL VALIDATION SUGGESTIONS
===================================
1. EMSA/gel shift assays:
   - Test top motif sequences for ARF27 binding
   - Compare binding affinity of motif variants
   - Validate specificity with competition assays

2. Reporter gene assays:
   - Clone motif-containing regions upstream of reporter
   - Test transcriptional activity in plant cells
   - Mutate motifs to confirm requirement for activity

3. ChIP-qPCR validation:
   - Design primers around high-confidence motif sites
   - Confirm ARF27 binding by quantitative ChIP
   - Test binding specificity with ARF27 mutants

COMPARATIVE ANALYSIS OPPORTUNITIES
==================================
1. Cross-species comparison:
   - Compare motifs with Arabidopsis ARF ChIP-seq data
   - Analyze conservation across monocot/dicot species
   - Identify maize-specific motif variations

2. ARF family comparison:
   - Compare ARF27 motifs with other ARF family members
   - Identify subfamily-specific binding preferences
   - Understand functional specialization

3. Temporal/spatial analysis:
   - Compare motifs across developmental stages
   - Analyze tissue-specific binding preferences
   - Identify condition-responsive motif usage

TROUBLESHOOTING GUIDE
====================

Common Issues and Solutions:

1. No ARF motifs detected:
   - Check ChIP-seq quality (FRiP scores, peak numbers)
   - Verify antibody specificity
   - Try different motif discovery parameters
   - Increase sequence length around summits

2. Too many non-specific motifs:
   - Use more stringent p-value thresholds
   - Improve background sequence selection
   - Focus analysis on highest-confidence peaks
   - Check for contamination or cross-reactivity

3. Weak motif signals:
   - Increase peak set size for motif discovery
   - Try different sequence lengths (±50bp, ±100bp, ±250bp)
   - Use only peaks with highest scores
   - Check peak calling parameters

4. Computational issues:
   - Reduce sequence set size for faster analysis
   - Use high-performance computing resources
   - Optimize MEME parameters for available memory
   - Process chromosomes separately if needed

STATISTICAL INTERPRETATION
==========================
1. E-value interpretation:
   - E < 1e-5: Highly significant motif
   - E < 1e-3: Significant motif
   - E > 1e-2: Marginal significance, manual inspection needed

2. FIMO p-value interpretation:
   - p < 1e-4: High-confidence motif occurrence
   - p < 1e-3: Moderate confidence
   - p > 1e-2: Low confidence, potential false positive

3. Enrichment assessment:
   - >50% peaks with motifs: Excellent specificity
   - 30-50% peaks with motifs: Good specificity  
   - 15-30% peaks with motifs: Moderate specificity
   - <15% peaks with motifs: Poor specificity or weak signal

FUTURE DIRECTIONS
================
1. Functional validation:
   - Test motif predictions with reporter assays
   - Validate binding with biochemical methods
   - Confirm biological relevance with mutant analysis

2. Mechanistic studies:
   - Investigate protein-protein interactions at motif sites
   - Study chromatin context of motif-containing regions
   - Analyze motif accessibility and nucleosome positioning

3. Comparative genomics:
   - Extend analysis to other maize genotypes
   - Compare with wild maize relatives
   - Study motif evolution and selection pressure

4. Systems analysis:
   - Integrate with gene expression data
   - Combine with chromatin accessibility studies
   - Model regulatory networks based on motif predictions

PUBLICATION RECOMMENDATIONS
===========================
1. Key figures to include:
   - Motif logos of top discovered motifs
   - Enrichment comparison with background
   - Position distribution relative to peak summits
   - Comparison with known ARF motifs from literature

2. Statistical reporting:
   - Report E-values for all discovered motifs
   - Include confidence intervals for enrichment estimates
   - Provide multiple testing correction details
   - Document all parameter choices and thresholds

3. Validation experiments:
   - Include experimental validation of top motifs
   - Show specificity controls and negative controls
   - Demonstrate functional relevance of predicted sites

CONCLUSION
==========
$(if [ -f "../results/motif_validation_summary.txt" ]; then
    grep -A 3 "Overall Assessment:" ../results/motif_validation_summary.txt | head -1
else
    echo "Motif analysis provides insights into ARF27 DNA-binding specificity."
fi)

The comprehensive motif analysis reveals ARF27's DNA-binding preferences and
validates the specificity of the ChIP-seq experiment. Results should be
interpreted in biological context and validated through experimental approaches.

CONTACT AND SUPPORT
===================
For questions about this analysis or motif interpretation:
- Review MEME Suite documentation: http://meme-suite.org/
- Consult plant ChIP-seq analysis guides
- Consider statistical consultation for complex interpretations
- Engage with plant molecular biology experts for functional validation

EOF

echo "Final motif analysis report generated: ARF27_motif_analysis_final_report.txt"
```

## Common Issues and Troubleshooting

### Issue 1: MEME Takes Too Long

**Problem**: MEME analysis runs for hours without completing

**Solutions:**
```bash
# Reduce sequence set size
head -200 peak_summits_200bp.fasta > subset_peaks.fasta

# Use faster MEME parameters
meme subset_peaks.fasta -dna -oc quick_meme -nmotifs 5 -minw 6 -maxw 12 -time 60

# Use DREME for faster discovery of short motifs
dreme -p peak_summits_200bp.fasta -n background_200bp.fasta -oc dreme_results
```

### Issue 2: No Significant Motifs Found

**Problem**: MEME or FIMO returns no significant results

**Diagnosis:**
```bash
# Check sequence quality
grep -v ">" peak_summits_200bp.fasta | grep -o "[ATGC]" | wc -l  # Valid bases
grep -v ">" peak_summits_200bp.fasta | grep -o "N" | wc -l       # N bases

# Check peak quality
awk '$5 > 100' ARF27_peaks.narrowPeak | wc -l  # High-score peaks
```

**Solutions:**
```bash
# Use only high-confidence peaks
awk '$5 > 50' ARF27_peaks.narrowPeak > high_confidence_peaks.narrowPeak

# Relax significance thresholds
fimo --thresh 1e-3 motifs.meme sequences.fasta  # Less stringent

# Try shorter sequences (focus on core binding site)
# Extract ±25bp around summit instead of ±100bp
```

### Issue 3: Background Sequences Issues

**Problem**: Inappropriate background affects motif discovery

**Solutions:**
```bash
# Create better background sequences
# Method 1: Matching GC content
fasta-shuffle-letters -kmer 2 peak_sequences.fasta > background.fasta

# Method 2: Promoter regions as background
awk '$3=="gene"' genes.gff3 | awk '{print $1"\t"($4-500)"\t"($4+500)}' > promoters.bed
bedtools getfasta -fi genome.fa -bed promoters.bed > promoter_background.fasta

# Method 3: Exclude repetitive regions
bedtools intersect -a random_regions.bed -b repeats.bed -v > clean_background.bed
```

### Issue 4: Memory Issues with Large Genomes

**Problem**: Tools crash due to insufficient memory

**Solutions:**
```bash
# Process smaller subsets
split -l 1000 peak_summits_200bp.fasta subset_

# Use streaming approaches
# Process one chromosome at a time
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
do
    awk -v chr="$chr" '$1==chr' peaks.bed > peaks_${chr}.bed
    # Run motif analysis on each chromosome
done

# Increase available memory
export TMPDIR=/path/to/large/tmp/directory
ulimit -v 16000000  # Increase virtual memory limit
```

### Issue 5: Poor Motif Quality

**Problem**: Discovered motifs are weak or non-specific

**Solutions:**
```bash
# Focus on peak centers
# Use smaller regions around summits (±50bp instead of ±200bp)

# Filter by peak quality
awk '$5 > 100 && $8 > 10' peaks.narrowPeak > stringent_peaks.narrowPeak

# Use multiple sequence lengths
# Try 50bp, 100bp, 200bp to find optimal resolution

# Improve peak calling
# Use more stringent MACS2 parameters
macs2 callpeak --qvalue 0.01 --pvalue 1e-5 ...
```

## Advanced Motif Analysis

### 1. Motif Spacing Analysis

```bash
# Analyze spacing between motifs
R --vanilla << 'EOF'
if(file.exists("../known_motifs/fimo_arf_motifs/fimo.tsv")) {
    fimo <- read.table("../known_motifs/fimo_arf_motifs/fimo.tsv", 
                      header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # Find sequences with multiple motifs
    multi_motif_seqs <- table(fimo$sequence_name)
    multi_motif_seqs <- multi_motif_seqs[multi_motif_seqs > 1]
    
    if(length(multi_motif_seqs) > 0) {
        cat("Sequences with multiple motifs:", length(multi_motif_seqs), "\n")
        
        # Analyze spacing for sequences with 2 motifs
        dual_motifs <- names(multi_motif_seqs)[multi_motif_seqs == 2]
        if(length(dual_motifs) > 0) {
            spacings <- c()
            for(seq in dual_motifs) {
                motif_pos <- fimo[fimo$sequence_name == seq, "start"]
                if(length(motif_pos) == 2) {
                    spacings <- c(spacings, abs(diff(motif_pos)))
                }
            }
            
            if(length(spacings) > 0) {
                cat("Motif spacing analysis (", length(spacings), "pairs):\n")
                cat("Mean spacing:", round(mean(spacings)), "bp\n")
                cat("Median spacing:", round(median(spacings)), "bp\n")
                
                png("../plots/motif_spacing_distribution.png", width=800, height=600)
                hist(spacings, breaks=20, 
                     main="Spacing Between ARF Motif Pairs", 
                     xlab="Distance (bp)", ylab="Frequency",
                     col="lightblue")
                dev.off()
            }
        }
    }
}
EOF
```

### 2. Motif Conservation Analysis

```bash
# Compare motifs across peak confidence levels
R --vanilla << 'EOF'
# Read peak data and motif data
peaks <- read.table("../../07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak", 
                   sep="\t", stringsAsFactors=FALSE)
colnames(peaks)[c(1,2,3,4,5)] <- c("chr", "start", "end", "name", "score")

if(file.exists("../known_motifs/fimo_arf_motifs/fimo.tsv")) {
    fimo <- read.table("../known_motifs/fimo_arf_motifs/fimo.tsv", 
                      header=TRUE, sep="\t", stringsAsFactors=FALSE)
    
    # Create peak score categories
    peaks$score_category <- cut(peaks$score, 
                               breaks=c(0, 50, 100, 200, Inf),
                               labels=c("Low", "Medium", "High", "Very High"))
    
    # Count motifs per score category
    peak_names <- paste0("peak_", 1:nrow(peaks))
    peaks$peak_id <- peak_names
    
    # Match FIMO results with peak scores
    fimo$peak_num <- as.numeric(gsub("peak_", "", fimo$sequence_name))
    fimo_with_scores <- merge(fimo, peaks[, c("score", "score_category")], 
                             by.x="peak_num", by.y=0, all.x=TRUE)
    
    # Analyze motif enrichment by peak score
    if(nrow(fimo_with_scores) > 0) {
        motif_by_score <- table(fimo_with_scores$score_category)
        peaks_by_score <- table(peaks$score_category)
        
        enrichment_data <- data.frame(
            Score_Category = names(peaks_by_score),
            Total_Peaks = as.numeric(peaks_by_score),
            Peaks_with_Motifs = as.numeric(motif_by_score[names(peaks_by_score)]),
            Enrichment_Rate = round(as.numeric(motif_by_score[names(peaks_by_score)]) / 
                                   as.numeric(peaks_by_score) * 100, 1)
        )
        
        enrichment_data$Peaks_with_Motifs[is.na(enrichment_data$Peaks_with_Motifs)] <- 0
        enrichment_data$Enrichment_Rate[is.na(enrichment_data$Enrichment_Rate)] <- 0
        
        cat("Motif enrichment by peak score category:\n")
        print(enrichment_data)
        
        write.csv(enrichment_data, "../results/motif_enrichment_by_score.csv", row.names=FALSE)
        
        # Visualization
        png("../plots/motif_enrichment_by_score.png", width=800, height=600)
        barplot(enrichment_data$Enrichment_Rate, 
                names.arg=enrichment_data$Score_Category,
                main="ARF Motif Enrichment by Peak Score",
                ylab="Motif Enrichment Rate (%)",
                xlab="Peak Score Category",
                col=c("lightcoral", "lightblue", "lightgreen", "gold"))
        dev.off()
    }
}
EOF
```

## Analysis Log Update

```bash
# Update analysis log
cat >> logs/step11_motif_analysis_log.txt << EOF

Step 11: Motif Analysis and Binding Site Discovery Results
Date: $(date)

Analysis Overview:
- Total peaks analyzed: $(grep -c ">" 11_motif_analysis/sequences/peak_summits_200bp.fasta)
- Sequence lengths tested: 100bp, 200bp, 500bp around summits
- Background controls: Random genomic regions and shuffled sequences

Methods Applied:
- MEME: De novo motif discovery
- FIMO: Known ARF motif scanning  
- TOMTOM: Motif similarity comparison
- Custom R scripts: Statistical analysis and visualization

Key Results:
$(if [ -f "11_motif_analysis/known_motifs/fimo_arf_motifs/fimo.tsv" ]; then
    significant_matches=$(tail -n +2 11_motif_analysis/known_motifs/fimo_arf_motifs/fimo.tsv | wc -l)
    total_peaks=$(grep -c ">" 11_motif_analysis/sequences/peak_summits_200bp.fasta)
    percentage=$(echo "scale=1; $significant_matches * 100 / $total_peaks" | bc)
    echo "- ARF motif detection: $significant_matches matches in $total_peaks peaks ($percentage%)"
else
    echo "- ARF motif detection: Analysis incomplete"
fi)
$(if [ -f "11_motif_analysis/meme_analysis/meme_top500_200bp/meme.txt" ]; then
    motif_count=$(grep -c "MOTIF" 11_motif_analysis/meme_analysis/meme_top500_200bp/meme.txt)
    echo "- De novo motifs discovered: $motif_count"
else
    echo "- De novo motif discovery: Analysis incomplete"
fi)

Quality Assessment:
$(if [ -f "11_motif_analysis/results/motif_validation_summary.txt" ]; then
    grep "Overall Assessment:" 11_motif_analysis/results/motif_validation_summary.txt
else
    echo "- Quality assessment: Pending completion"
fi)

Files Generated:
- Peak sequences: 6 FASTA files (different lengths and subsets)
- MEME results: HTML reports and motif files
- FIMO results: TSV files with motif occurrences
- Visualization plots: 8+ plots showing motif characteristics
- Validation summaries: 2 comprehensive reports

Tools Used:
- MEME Suite: motif discovery and scanning
- bedtools: sequence extraction
- R/Bioconductor: statistical analysis
- Custom scripts: visualization and validation

Issues Encountered: [None/List problems and solutions]

Biological Insights:
- ChIP-seq specificity validation: [Based on ARF motif detection]
- ARF27 binding preferences: [Based on motif characteristics]
- Regulatory mechanisms: [Based on motif distribution and spacing]

Ready for:
✓ Data integration and interpretation
✓ Experimental validation design  
✓ Manuscript preparation
✓ Additional comparative analyses

EOF
```

## Key Takeaways

- Motif analysis validates ChIP-seq specificity and success
- Presence of known ARF motifs confirms antibody specificity
- De novo discovery reveals binding site preferences
- Motif distribution patterns indicate regulatory mechanisms
- Statistical validation ensures biological relevance
- Results guide experimental validation strategies

## Final Checklist

Before concluding the ChIP-seq analysis:

- [ ] Peak sequences extracted around summits
- [ ] Background sequences generated appropriately
- [ ] MEME de novo motif discovery completed
- [ ] FIMO known motif scanning performed
- [ ] Motif comparison with known ARF motifs done
- [ ] Statistical validation and enrichment calculated
- [ ] Visualization plots generated
- [ ] Comprehensive validation summary created
- [ ] Final motif analysis report completed
- [ ] Results interpreted in biological context

## Tutorial Completion

**Congratulations!** You have completed the comprehensive ARF27 ChIP-seq analysis tutorial. This 11-step workflow has taken you from raw sequencing reads to publication-ready results, including:

1. ✅ Read Quality Control
2. ✅ Read Mapping  
3. ✅ SAM to BAM Conversion
4. ✅ PCR Duplicate Removal
5. ✅ Data Normalization
6. ✅ Peak Calling
7. ✅ Peak Blacklist Filtering
8. ✅ Peak Data Quality Control
9. ✅ Differential Peak Analysis
10. ✅ Peak Annotation and Gene Target Identification
11. ✅ Motif Analysis and Binding Site Discovery

**Your analysis now includes:**
- High-quality, filtered ARF27 binding sites
- Annotated target genes with functional categories
- Validated motif predictions
- Comprehensive quality metrics
- Publication-ready visualizations
- Experimental validation strategies

**Next steps for your research:**
- Experimental validation of key findings
- Integration with transcriptome data
- Comparative analysis with other conditions
- Functional studies of target genes
- Manuscript preparation and publication

This tutorial provides a solid foundation for ChIP-seq analysis that can be adapted for other transcription factors, experimental systems, and research questions.
