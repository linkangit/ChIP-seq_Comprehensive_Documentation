# Step 10: Peak Annotation and Gene Target Identification

## Objective
Annotate ARF27 binding peaks relative to genes and genomic features to identify potential target genes, understand the regulatory landscape, and gain biological insights into ARF27 function in maize.

## Why Peak Annotation is Critical

### Biological Questions Addressed:
- **Target gene identification**: Which genes are directly regulated by ARF27?
- **Regulatory mechanisms**: Does ARF27 bind promoters, enhancers, or gene bodies?
- **Functional categories**: What biological processes are controlled by ARF27?
- **Evolutionary conservation**: Are targets conserved across species?

### Types of Annotation:
1. **Distance-based**: Assign peaks to nearest genes
2. **Feature-based**: Classify peaks by genomic regions (promoter, exon, intron, intergenic)
3. **Functional**: Gene ontology and pathway enrichment
4. **Comparative**: Cross-species target conservation

### ARF27-Specific Expectations:
- **Promoter binding**: Expected for transcriptional regulation
- **Auxin-responsive genes**: Targets should be enriched for auxin signaling
- **Developmental genes**: ARF27 regulates plant development
- **Stress response genes**: ARFs respond to environmental conditions

## Step-by-Step Peak Annotation

### 1. Set Up Annotation Directory

```bash
# Navigate to analysis directory
cd chipseq_arf27_analysis

# Create annotation directory structure
mkdir -p 10_peak_annotation
mkdir -p 10_peak_annotation/gene_annotation
mkdir -p 10_peak_annotation/feature_annotation
mkdir -p 10_peak_annotation/functional_annotation
mkdir -p 10_peak_annotation/target_genes
mkdir -p 10_peak_annotation/plots
mkdir -p 10_peak_annotation/databases
```

### 2. Download and Prepare Gene Annotations

Download maize gene annotations and prepare for analysis:

```bash
# Download maize B73-v4 gene annotation
cd 10_peak_annotation/databases

# Option 1: Download from MaizeGDB
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3

# Option 2: Download from Ensembl Plants
# wget http://ftp.ensemblgenomes.org/pub/plants/release-52/gff3/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gff3.gz
# gunzip Zea_mays.Zm-B73-REFERENCE-NAM-5.0.52.gff3.gz

# Rename for clarity
mv Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3 maize_B73_v4_genes.gff3

# Extract different feature types
echo "Extracting genomic features..."

# Extract gene features
awk '$3=="gene"' maize_B73_v4_genes.gff3 > ../gene_annotation/maize_genes.gff3

# Extract mRNA/transcript features
awk '$3=="mRNA" || $3=="transcript"' maize_B73_v4_genes.gff3 > ../gene_annotation/maize_transcripts.gff3

# Extract exon features
awk '$3=="exon"' maize_B73_v4_genes.gff3 > ../gene_annotation/maize_exons.gff3

# Convert to BED format for easier processing
echo "Converting to BED format..."

# Genes to BED
awk '$3=="gene" {
    match($9, /ID=([^;]+)/, id);
    match($9, /Name=([^;]+)/, name);
    print $1"\t"($4-1)"\t"$5"\t"id[1]"\t"$6"\t"$7"\t"(name[1] ? name[1] : id[1])
}' maize_B73_v4_genes.gff3 > ../gene_annotation/maize_genes.bed

# Create promoter regions (2kb upstream of TSS)
awk '$3=="gene" {
    match($9, /ID=([^;]+)/, id);
    if($7=="+") {
        start = ($4-2000 > 0) ? $4-2000 : 1;
        end = $4;
    } else {
        start = $5;
        end = $5+2000;
    }
    print $1"\t"start"\t"end"\t"id[1]"_promoter\t"$6"\t"$7
}' maize_B73_v4_genes.gff3 > ../gene_annotation/maize_promoters.bed

echo "Gene annotation files prepared."
```

### 3. Install Peak Annotation Tools

```bash
# Install R packages for peak annotation
R --vanilla << 'EOF'
# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "ChIPseeker",
    "TxDb.Zmays.NCBI.B73",  # If available
    "GenomicFeatures",
    "rtracklayer",
    "GenomicRanges",
    "clusterProfiler",
    "org.Zm.eg.db",  # If available
    "DOSE",
    "enrichplot"
))

# Install CRAN packages
install.packages(c(
    "ggplot2",
    "dplyr",
    "VennDiagram",
    "pheatmap"
))

# Check installations
library(ChIPseeker)
library(GenomicFeatures)
sessionInfo()
EOF
```

### 4. Distance-Based Gene Annotation

Assign peaks to their nearest genes:

```bash
# Annotate peaks with nearest genes
echo "Performing distance-based gene annotation..."

cd ../../  # Return to main analysis directory

# Use bedtools to find nearest genes
bedtools closest \
    -a 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak \
    -b 10_peak_annotation/gene_annotation/maize_genes.bed \
    -d > 10_peak_annotation/gene_annotation/peaks_nearest_genes.txt

# Process results to create clean annotation
R --vanilla << 'EOF'
# Read nearest gene annotation
nearest_genes <- read.table("10_peak_annotation/gene_annotation/peaks_nearest_genes.txt", 
                           sep="\t", stringsAsFactors=FALSE)

# Create column names
colnames(nearest_genes) <- c(
    "peak_chr", "peak_start", "peak_end", "peak_name", "peak_score", "peak_strand",
    "peak_signalValue", "peak_pValue", "peak_qValue", "peak_summit",
    "gene_chr", "gene_start", "gene_end", "gene_id", "gene_score", "gene_strand", "gene_name",
    "distance"
)

# Clean up annotation
annotation_clean <- data.frame(
    peak_id = nearest_genes$peak_name,
    peak_chr = nearest_genes$peak_chr,
    peak_start = nearest_genes$peak_start,
    peak_end = nearest_genes$peak_end,
    peak_score = nearest_genes$peak_score,
    gene_id = nearest_genes$gene_id,
    gene_name = nearest_genes$gene_name,
    distance_to_gene = nearest_genes$distance,
    peak_location = ifelse(nearest_genes$distance == 0, "gene_body",
                          ifelse(abs(nearest_genes$distance) <= 2000, "promoter",
                                ifelse(abs(nearest_genes$distance) <= 10000, "proximal",
                                      "distal"))),
    stringsAsFactors = FALSE
)

# Save annotation
write.csv(annotation_clean, 
          "10_peak_annotation/gene_annotation/peak_gene_annotation.csv", 
          row.names = FALSE)

# Generate summary statistics
cat("Peak Annotation Summary:\n")
cat("Total peaks annotated:", nrow(annotation_clean), "\n")
cat("Gene body peaks:", sum(annotation_clean$peak_location == "gene_body"), "\n")
cat("Promoter peaks (≤2kb):", sum(annotation_clean$peak_location == "promoter"), "\n")
cat("Proximal peaks (2-10kb):", sum(annotation_clean$peak_location == "proximal"), "\n")
cat("Distal peaks (>10kb):", sum(annotation_clean$peak_location == "distal"), "\n")

# Create distance distribution plot
png("10_peak_annotation/plots/distance_distribution.png", width=800, height=600)
hist(log10(abs(annotation_clean$distance_to_gene) + 1), breaks=50,
     main="Distance to Nearest Gene Distribution",
     xlab="Log10(Distance to TSS + 1)",
     ylab="Number of Peaks",
     col="lightblue")
abline(v=log10(2001), col="red", lty=2, lwd=2)
text(log10(2001), max(hist(log10(abs(annotation_clean$distance_to_gene) + 1), plot=FALSE)$counts)/2, 
     "2kb", pos=4, col="red")
dev.off()

# Create location pie chart
location_counts <- table(annotation_clean$peak_location)
png("10_peak_annotation/plots/peak_location_pie.png", width=600, height=600)
pie(location_counts, main="Peak Location Distribution", 
    col=c("lightcoral", "lightblue", "lightgreen", "lightyellow"))
dev.off()

cat("Distance-based annotation completed.\n")
EOF
```

### 5. Feature-Based Annotation Using ChIPseeker

Comprehensive genomic feature annotation:

```bash
# Use ChIPseeker for detailed annotation
R --vanilla << 'EOF'
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)
library(ggplot2)

# Create TxDb object from GFF3
# Note: This may take several minutes for large genomes
cat("Creating transcript database from GFF3...\n")

# Read GFF3 file
gff_file <- "10_peak_annotation/databases/maize_B73_v4_genes.gff3"

# Create TxDb object
txdb <- makeTxDbFromGFF(gff_file, format="gff3")

# Save TxDb for future use
saveDb(txdb, "10_peak_annotation/databases/maize_B73_v4_txdb.sqlite")

# Read peak file
peak_file <- "07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak"
peaks <- readPeakFile(peak_file)

# Perform ChIPseeker annotation
cat("Annotating peaks with ChIPseeker...\n")
peak_anno <- annotatePeak(peaks, 
                         tssRegion=c(-2000, 2000),
                         TxDb=txdb,
                         annoDb=NULL)  # No org.db available for maize

# Save annotation results
anno_df <- as.data.frame(peak_anno)
write.csv(anno_df, "10_peak_annotation/feature_annotation/chipseeker_annotation.csv", 
          row.names=FALSE)

# Create visualization plots
cat("Creating ChIPseeker visualization plots...\n")

# TSS plot
png("10_peak_annotation/plots/tss_plot.png", width=800, height=600)
plotAvgProf(peak_anno, xlim=c(-2000, 2000), xlab="Distance from TSS (bp)", 
            ylab="Read Count Frequency")
dev.off()

# Annotation pie chart
png("10_peak_annotation/plots/annotation_pie.png", width=800, height=600)
plotAnnoPie(peak_anno)
dev.off()

# Distance to TSS
png("10_peak_annotation/plots/distance_to_tss.png", width=800, height=600)
plotDistToTSS(peak_anno, title="Distance to TSS")
dev.off()

# Annotation bar plot
png("10_peak_annotation/plots/annotation_bar.png", width=1000, height=600)
plotAnnoBar(peak_anno)
dev.off()

# Create summary
annotation_summary <- data.frame(
    Feature = names(table(anno_df$annotation)),
    Count = as.numeric(table(anno_df$annotation)),
    Percentage = round(as.numeric(table(anno_df$annotation))/nrow(anno_df)*100, 1)
)

write.csv(annotation_summary, 
          "10_peak_annotation/feature_annotation/annotation_summary.csv", 
          row.names=FALSE)

cat("ChIPseeker annotation completed.\n")
print(annotation_summary)
EOF
```

### 6. Identify High-Confidence Target Genes

Extract the most likely direct targets:

```bash
# Identify high-confidence ARF27 target genes
echo "Identifying high-confidence target genes..."

R --vanilla << 'EOF'
# Read peak annotation
if(file.exists("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")) {
    annotation <- read.csv("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")
    
    # Define high-confidence targets based on multiple criteria
    high_conf_targets <- annotation[
        (annotation$peak_location %in% c("gene_body", "promoter")) &
        (annotation$peak_score > 50) &  # High-quality peaks
        (abs(annotation$distance_to_gene) <= 5000),  # Close to gene
    ]
    
    # Remove duplicate gene targets (keep highest scoring peak per gene)
    high_conf_targets <- high_conf_targets[order(high_conf_targets$gene_id, 
                                                -high_conf_targets$peak_score), ]
    high_conf_unique <- high_conf_targets[!duplicated(high_conf_targets$gene_id), ]
    
    # Create target gene list
    target_genes <- data.frame(
        gene_id = high_conf_unique$gene_id,
        gene_name = high_conf_unique$gene_name,
        peak_score = high_conf_unique$peak_score,
        distance_to_gene = high_conf_unique$distance_to_gene,
        peak_location = high_conf_unique$peak_location,
        stringsAsFactors = FALSE
    )
    
    # Sort by peak score (highest confidence first)
    target_genes <- target_genes[order(-target_genes$peak_score), ]
    
    # Save target gene list
    write.csv(target_genes, 
              "10_peak_annotation/target_genes/high_confidence_targets.csv", 
              row.names = FALSE)
    
    cat("High-confidence target identification:\n")
    cat("Total annotated peaks:", nrow(annotation), "\n")
    cat("High-confidence target genes:", nrow(target_genes), "\n")
    cat("Promoter targets:", sum(target_genes$peak_location == "promoter"), "\n")
    cat("Gene body targets:", sum(target_genes$peak_location == "gene_body"), "\n")
    
    # Create score distribution plot
    png("10_peak_annotation/plots/target_gene_scores.png", width=800, height=600)
    hist(target_genes$peak_score, breaks=30,
         main="ARF27 Target Gene Peak Scores",
         xlab="Peak Score",
         ylab="Number of Target Genes",
         col="lightgreen")
    dev.off()
    
    # Create distance distribution for targets
    png("10_peak_annotation/plots/target_distance_distribution.png", width=800, height=600)
    boxplot(abs(target_genes$distance_to_gene) ~ target_genes$peak_location,
            main="Distance Distribution by Peak Location",
            xlab="Peak Location", ylab="Distance to Gene (bp)",
            col=c("lightblue", "lightcoral"))
    dev.off()
    
} else {
    cat("Peak annotation file not found. Please run distance-based annotation first.\n")
}
EOF
```

### 7. Functional Enrichment Analysis

Analyze the biological functions of ARF27 target genes:

```bash
# Perform functional enrichment analysis
echo "Performing functional enrichment analysis..."

# First, we need to prepare gene lists and functional annotations
# Since maize-specific databases may be limited, we'll use plant-focused approaches

R --vanilla << 'EOF'
library(dplyr)

# Read target genes
if(file.exists("10_peak_annotation/target_genes/high_confidence_targets.csv")) {
    targets <- read.csv("10_peak_annotation/target_genes/high_confidence_targets.csv")
    
    # Extract gene IDs for enrichment analysis
    target_gene_ids <- targets$gene_id
    
    # Save gene lists for external analysis
    writeLines(target_gene_ids, "10_peak_annotation/functional_annotation/target_gene_ids.txt")
    
    # Create background gene list (all genes with peaks nearby)
    if(file.exists("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")) {
        all_annotation <- read.csv("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")
        background_genes <- unique(all_annotation$gene_id)
        writeLines(background_genes, "10_peak_annotation/functional_annotation/background_gene_ids.txt")
        
        cat("Gene lists created for functional analysis:\n")
        cat("Target genes:", length(target_gene_ids), "\n")
        cat("Background genes:", length(background_genes), "\n")
    }
    
    # Basic functional categories based on gene names/descriptions
    # This is a simplified approach - real analysis would use GO terms or pathway databases
    
    # Look for auxin-related genes (simplified pattern matching)
    auxin_related <- grep("auxin|AUX|IAA|ARF", targets$gene_name, ignore.case=TRUE)
    stress_related <- grep("stress|drought|salt|heat|cold", targets$gene_name, ignore.case=TRUE)
    development_related <- grep("development|growth|embryo|root|shoot", targets$gene_name, ignore.case=TRUE)
    
    functional_summary <- data.frame(
        Category = c("Total targets", "Auxin-related", "Stress-related", "Development-related"),
        Count = c(nrow(targets), length(auxin_related), length(stress_related), length(development_related)),
        Percentage = c(100, 
                      round(length(auxin_related)/nrow(targets)*100, 1),
                      round(length(stress_related)/nrow(targets)*100, 1),
                      round(length(development_related)/nrow(targets)*100, 1))
    )
    
    write.csv(functional_summary, 
              "10_peak_annotation/functional_annotation/functional_summary.csv", 
              row.names=FALSE)
    
    print(functional_summary)
    
    # Create functional category plot
    png("10_peak_annotation/plots/functional_categories.png", width=800, height=600)
    functional_counts <- functional_summary$Count[2:4]
    names(functional_counts) <- functional_summary$Category[2:4]
    barplot(functional_counts, 
            main="Functional Categories of ARF27 Target Genes",
            ylab="Number of Genes",
            col=c("lightcoral", "lightblue", "lightgreen"),
            las=2)
    dev.off()
    
} else {
    cat("Target gene file not found. Please run target identification first.\n")
}

# Instructions for external enrichment analysis
cat("\nFor comprehensive functional enrichment analysis:\n")
cat("1. Use target_gene_ids.txt with online tools like:\n")
cat("   - AgriGO (http://systemsbiology.cau.edu.cn/agriGOv2/)\n")
cat("   - PlantGSEA (http://structuralbiology.cau.edu.cn/PlantGSEA/)\n")
cat("   - MaizeGDB tools\n")
cat("2. Upload background_gene_ids.txt as background set\n")
cat("3. Use Gene Ontology or KEGG pathway analysis\n")
EOF
```

### 8. Compare with Known ARF Targets

Compare with literature-known ARF targets:

```bash
# Create comparison with known auxin-responsive genes
echo "Comparing with known auxin-responsive genes..."

# Create a list of known auxin-responsive genes (simplified example)
cat > 10_peak_annotation/functional_annotation/known_auxin_genes.txt << 'EOF'
Zm00001eb000010
Zm00001eb000050
Zm00001eb000100
Zm00001eb000150
Zm00001eb000200
EOF

R --vanilla << 'EOF'
# Read known auxin genes and target genes
known_auxin <- readLines("10_peak_annotation/functional_annotation/known_auxin_genes.txt")

if(file.exists("10_peak_annotation/target_genes/high_confidence_targets.csv")) {
    targets <- read.csv("10_peak_annotation/target_genes/high_confidence_targets.csv")
    target_ids <- targets$gene_id
    
    # Find overlap with known auxin genes
    overlap <- intersect(target_ids, known_auxin)
    
    cat("Comparison with known auxin-responsive genes:\n")
    cat("Known auxin genes:", length(known_auxin), "\n")
    cat("ARF27 targets:", length(target_ids), "\n")
    cat("Overlap:", length(overlap), "\n")
    cat("Overlap percentage:", round(length(overlap)/length(known_auxin)*100, 1), "%\n")
    
    if(length(overlap) > 0) {
        cat("Overlapping genes:\n")
        print(overlap)
        
        # Save overlapping genes
        writeLines(overlap, "10_peak_annotation/functional_annotation/overlap_with_known_auxin.txt")
    }
    
    # Create Venn diagram if overlap exists
    if(length(overlap) > 0) {
        library(VennDiagram)
        
        venn.plot <- draw.pairwise.venn(
            area1 = length(known_auxin),
            area2 = length(target_ids),
            cross.area = length(overlap),
            category = c("Known Auxin Genes", "ARF27 Targets"),
            fill = c("lightblue", "lightgreen"),
            alpha = 0.5,
            filename = NULL
        )
        
        png("10_peak_annotation/plots/known_targets_overlap.png", width=600, height=600)
        grid.draw(venn.plot)
        dev.off()
    }
}
EOF
```

### 9. Create Target Gene Browser Tracks

Prepare files for genome browser visualization:

```bash
# Create browser tracks for target genes
echo "Creating browser tracks for target genes..."

# Create BED file of target gene regions
R --vanilla << 'EOF'
if(file.exists("10_peak_annotation/target_genes/high_confidence_targets.csv")) {
    targets <- read.csv("10_peak_annotation/target_genes/high_confidence_targets.csv")
    
    # Read full gene annotation to get coordinates
    if(file.exists("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")) {
        annotation <- read.csv("10_peak_annotation/gene_annotation/peak_gene_annotation.csv")
        
        # Get target gene coordinates
        target_coords <- annotation[annotation$gene_id %in% targets$gene_id, 
                                   c("peak_chr", "peak_start", "peak_end", "gene_id")]
        
        # Remove duplicates and sort
        target_coords <- target_coords[!duplicated(target_coords$gene_id), ]
        target_coords <- target_coords[order(target_coords$peak_chr, target_coords$peak_start), ]
        
        # Create BED file
        bed_format <- data.frame(
            chr = target_coords$peak_chr,
            start = target_coords$peak_start,
            end = target_coords$peak_end,
            name = target_coords$gene_id,
            score = 1000,
            strand = "."
        )
        
        write.table(bed_format, 
                   "10_peak_annotation/target_genes/target_genes.bed",
                   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
        
        cat("Browser track created: target_genes.bed\n")
    }
}
EOF

# Create GTF file for target genes (if needed)
awk 'BEGIN{OFS="\t"} FNR==NR{targets[$1]=1; next} 
     $3=="gene" && ($9 ~ /ID=([^;]+)/) {
         match($9, /ID=([^;]+)/, id);
         if(id[1] in targets) print
     }' 10_peak_annotation/target_genes/high_confidence_targets.csv \
       10_peak_annotation/databases/maize_B73_v4_genes.gff3 \
       > 10_peak_annotation/target_genes/target_genes.gtf

echo "Target gene tracks prepared for genome browser visualization."
```

### 10. Generate Comprehensive Annotation Report

Create a final comprehensive report:

```bash
# Generate comprehensive annotation report
echo "Generating comprehensive annotation report..."

cat > 10_peak_annotation/ARF27_peak_annotation_report.txt << EOF
ARF27 ChIP-seq Peak Annotation Report
====================================
Analysis Date: $(date)
Genome Assembly: Maize B73-v4

ANNOTATION OVERVIEW
==================
Total peaks analyzed: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
Successfully annotated peaks: $(if [ -f "10_peak_annotation/gene_annotation/peak_gene_annotation.csv" ]; then tail -n +2 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l; else echo "N/A"; fi)

PEAK LOCATION DISTRIBUTION
=========================
$(if [ -f "10_peak_annotation/feature_annotation/annotation_summary.csv" ]; then
    echo "ChIPseeker annotation:"
    cat 10_peak_annotation/feature_annotation/annotation_summary.csv
else
    echo "ChIPseeker annotation not available"
fi)

DISTANCE-BASED ANNOTATION
========================
$(if [ -f "10_peak_annotation/gene_annotation/peak_gene_annotation.csv" ]; then
    echo "Gene body peaks: $(awk -F',' '$9=="gene_body"' 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)"
    echo "Promoter peaks (≤2kb): $(awk -F',' '$9=="promoter"' 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)"
    echo "Proximal peaks (2-10kb): $(awk -F',' '$9=="proximal"' 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)"
    echo "Distal peaks (>10kb): $(awk -F',' '$9=="distal"' 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)"
else
    echo "Distance-based annotation not available"
fi)

TARGET GENE IDENTIFICATION
=========================
$(if [ -f "10_peak_annotation/target_genes/high_confidence_targets.csv" ]; then
    echo "High-confidence target genes: $(tail -n +2 10_peak_annotation/target_genes/high_confidence_targets.csv | wc -l)"
    echo "Promoter targets: $(awk -F',' '$5=="promoter"' 10_peak_annotation/target_genes/high_confidence_targets.csv | wc -l)"
    echo "Gene body targets: $(awk -F',' '$5=="gene_body"' 10_peak_annotation/target_genes/high_confidence_targets.csv | wc -l)"
else
    echo "Target gene identification not completed"
fi)

FUNCTIONAL ANALYSIS
==================
$(if [ -f "10_peak_annotation/functional_annotation/functional_summary.csv" ]; then
    cat 10_peak_annotation/functional_annotation/functional_summary.csv
else
    echo "Functional analysis not completed"
fi)

FILES GENERATED
===============
1. Gene annotation files: gene_annotation/
2. Feature annotation files: feature_annotation/
3. Target gene lists: target_genes/
4. Functional analysis: functional_annotation/
5. Visualization plots: plots/
6. Browser tracks: target_genes/*.bed, *.gtf

QUALITY METRICS
===============
Annotation success rate: $(if [ -f "10_peak_annotation/gene_annotation/peak_gene_annotation.csv" ]; then
    total_peaks=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
    annotated_peaks=$(tail -n +2 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)
    echo "scale=1; $annotated_peaks * 100 / $total_peaks" | bc
else
    echo "N/A"
fi)%

Promoter enrichment: $(if [ -f "10_peak_annotation/gene_annotation/peak_gene_annotation.csv" ]; then
    total=$(tail -n +2 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)
    promoter=$(awk -F',' '$9=="promoter"' 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l)
    echo "scale=1; $promoter * 100 / $total" | bc
else
    echo "N/A"
fi)%

BIOLOGICAL INSIGHTS
==================
1. ARF27 binding pattern: [Interpret based on location distribution]
2. Target gene categories: [Based on functional analysis]
3. Regulatory mechanisms: [Promoter vs gene body binding]
4. Comparison with known targets: [If overlap analysis completed]

NEXT STEPS
==========
✓ Peak annotation completed
→ Ready for Step 11: Motif Analysis
→ Use target gene lists for:
  - External functional enrichment analysis
  - qRT-PCR validation
  - Comparative analysis with other ARF family members

EXTERNAL ANALYSIS RECOMMENDATIONS
===============================
1. Upload gene lists to AgriGO or PlantGSEA for GO enrichment
2. Check MaizeGDB for detailed gene information
3. Compare with published auxin-responsive gene sets
4. Validate key targets with qRT-PCR

EOF

echo "Comprehensive annotation report generated: 10_peak_annotation/ARF27_peak_annotation_report.txt"
```

## Analysis Log Update

```bash
# Update analysis log
cat >> logs/step10_peak_annotation_log.txt << EOF

Step 10: Peak Annotation and Gene Target Identification Results
Date: $(date)

Annotation Overview:
- Total peaks annotated: $(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
- Annotation success rate: $(if [ -f "10_peak_annotation/gene_annotation/peak_gene_annotation.csv" ]; then
    total_peaks=$(wc -l < 07_blacklists/filtered_peaks/ARF27_ChIP_merged_peaks_filtered.narrowPeak)
    annotated_peaks=$(tail -n +2 10_peak_annotation/gene_annotation/peak_gene_annotation.csv | wc -l
