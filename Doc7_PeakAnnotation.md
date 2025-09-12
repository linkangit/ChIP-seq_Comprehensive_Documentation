# Document 7: Peak Annotation and Functional Analysis

## Understanding Peak Annotation

Peak annotation is the process of determining where your ARF27 binding sites are located relative to genes and other genomic features. This step transforms a list of genomic coordinates into biological insights about which genes ARF27 regulates and what cellular processes it controls.

Think of annotation as adding biological context to your data. Instead of just knowing that ARF27 binds at position "chr1:1,234,567-1,235,000", you'll learn that this peak is:
- In the promoter of gene Zm00001d002345
- 500 bp upstream of the transcription start site
- Associated with a gene involved in auxin response

### The Biology of Transcription Factor Binding Locations

Where a transcription factor binds relative to genes determines how it functions:

**Promoter Binding** (within 2kb of transcription start sites):
- Most direct mechanism of gene regulation
- Can activate or repress transcription
- Often represents primary regulatory targets

**Gene Body Binding** (within exons or introns):
- May affect transcription elongation
- Could influence alternative splicing
- Sometimes represents indirect or cooperative binding

**Intergenic Binding** (between genes):
- May represent enhancer or silencer elements
- Could regulate distant genes through chromatin looping
- Might control non-coding RNAs

**Enhancer Binding** (distal regulatory elements):
- Can regulate genes from great distances
- Often shows tissue-specific activity
- Important for fine-tuning gene expression

### Why Functional Analysis Matters

Simply knowing where ARF27 binds isn't enough - we need to understand what biological processes it controls. Functional analysis reveals:
- Which cellular pathways are regulated by ARF27
- How ARF27 fits into auxin signaling networks
- What developmental processes depend on ARF27 function
- Which genes are primary vs. secondary targets

## Setting Up R for Annotation Analysis

We'll use R and Bioconductor packages for comprehensive annotation and functional analysis:

```bash
# Make sure we're in the right directory
cd chipseq_analysis

# Start R for annotation analysis
R
```

## Step 1: Load Required Libraries

```r
# Load required libraries for peak annotation
library(ChIPseeker)
library(GenomicFeatures)
library(rtracklayer)
library(org.Zm.eg.db)
library(clusterProfiler)
library(ggplot2)
library(tidyverse)

# Set up output directories
dir.create("annotation", showWarnings = FALSE)
dir.create("functional_analysis", showWarnings = FALSE)
```

### Understanding the R Packages

**ChIPseeker**: Specialized for ChIP-seq peak annotation
- Automatically annotates peaks relative to genes
- Creates publication-ready visualizations
- Handles complex genomic feature overlaps

**GenomicFeatures**: Infrastructure for genomic annotations
- Converts GTF files to database format
- Defines promoter and gene body regions
- Provides efficient genomic coordinate operations

**clusterProfiler**: Modern functional enrichment analysis
- Gene Ontology and pathway analysis
- Publication-ready plots
- Handles multiple testing correction properly

**org.Zm.eg.db**: Maize-specific annotation database
- Gene symbols, descriptions, and functional annotations
- Gene Ontology mappings for maize genes
- Cross-references between different gene identifiers

## Step 2: Create Transcript Database

```r
# Create transcript database from the GTF file
cat("Building transcript database from GTF file...\n")
txdb <- makeTxDbFromGFF("reference_genome/maize_annotations.gtf", 
                        format = "gtf",
                        organism = "Zea mays")

# Save the database for future use
saveDb(txdb, "reference_genome/maize_txdb.sqlite")

# Define promoter regions (2kb upstream, 1kb downstream of TSS)
promoter_regions <- getPromoters(TxDb = txdb, 
                                 upstream = 2000, 
                                 downstream = 1000)

cat("Transcript database created successfully!\n")
```

### Understanding Promoter Definition

**Why 2kb upstream, 1kb downstream?**
- Most transcription factor binding sites occur within 2kb of transcription start sites
- Some binding sites occur just downstream of TSS (within transcribed region)
- This definition captures most promoter-associated regulatory elements
- Can be adjusted based on your research questions

## Step 3: Annotate Peaks for Each Sample

### Annotate First Sample

```r
# Read and annotate peaks for sample 1
cat("Annotating peaks for sample1_ChIP...\n")

# Read peak file
sample1_peaks <- readPeakFile("peaks/sample1_ChIP_peaks.narrowPeak")

# Perform annotation
sample1_annotation <- annotatePeak(sample1_peaks, 
                                   TxDb = txdb,
                                   tssRegion = c(-2000, 1000),
                                   verbose = FALSE)

# Save detailed annotation results
sample1_df <- as.data.frame(sample1_annotation)
write.table(sample1_df, 
            "annotation/sample1_ChIP_detailed_annotation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("Sample1 annotation completed!\n")
```

### Annotate Additional Samples

```r
# Annotate sample 2
cat("Annotating peaks for sample2_ChIP...\n")
sample2_peaks <- readPeakFile("peaks/sample2_ChIP_peaks.narrowPeak")
sample2_annotation <- annotatePeak(sample2_peaks, 
                                   TxDb = txdb,
                                   tssRegion = c(-2000, 1000),
                                   verbose = FALSE)

sample2_df <- as.data.frame(sample2_annotation)
write.table(sample2_df, 
            "annotation/sample2_ChIP_detailed_annotation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

# Annotate sample 3
cat("Annotating peaks for sample3_ChIP...\n")
sample3_peaks <- readPeakFile("peaks/sample3_ChIP_peaks.narrowPeak")
sample3_annotation <- annotatePeak(sample3_peaks, 
                                   TxDb = txdb,
                                   tssRegion = c(-2000, 1000),
                                   verbose = FALSE)

sample3_df <- as.data.frame(sample3_annotation)
write.table(sample3_df, 
            "annotation/sample3_ChIP_detailed_annotation.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("All samples annotated!\n")
```

## Step 4: Create Annotation Visualizations

### Individual Sample Plots

```r
# Create annotation plots for sample 1
pdf("annotation/sample1_ChIP_annotation_plots.pdf", width = 10, height = 8)

# Peak distribution pie chart
plotAnnoPie(sample1_annotation, main = "Sample1 Peak Distribution")

# Distance to TSS distribution
plotDistToTSS(sample1_annotation, 
              title = "Distance to TSS - Sample1",
              xlab = "Distance to TSS (bp)")

# Annotation bar plot
plotAnnoBar(sample1_annotation, title = "Sample1 Peak Annotation")

dev.off()

cat("Sample1 plots created!\n")
```

### Compare Multiple Samples

```r
# Create comparison plots if you have multiple samples
if (exists("sample2_annotation") && exists("sample3_annotation")) {
    
    # Combine annotations for comparison
    peak_list <- list(
        Sample1 = sample1_annotation,
        Sample2 = sample2_annotation,
        Sample3 = sample3_annotation
    )
    
    # Create comparison plots
    pdf("annotation/multi_sample_comparison.pdf", width = 12, height = 8)
    
    # Compare peak distributions
    plotAnnoBar(peak_list, title = "Peak Annotation Comparison")
    
    # Compare distance to TSS
    plotDistToTSS(peak_list, title = "Distance to TSS Comparison")
    
    dev.off()
    
    cat("Comparison plots created!\n")
}
```

## Step 5: Create Annotation Summary

```r
# Create summary statistics for each sample
create_annotation_summary <- function(annotation_df, sample_name) {
    
    # Count annotations by type
    annotation_counts <- table(annotation_df$annotation)
    
    # Calculate percentages
    total_peaks <- nrow(annotation_df)
    annotation_summary <- data.frame(
        Annotation = names(annotation_counts),
        Count = as.numeric(annotation_counts),
        Percentage = round(as.numeric(annotation_counts) / total_peaks * 100, 1)
    )
    
    # Add sample information
    annotation_summary$Sample <- sample_name
    annotation_summary$Total_Peaks <- total_peaks
    
    return(annotation_summary)
}

# Generate summaries for each sample
sample1_summary <- create_annotation_summary(sample1_df, "Sample1")
write.table(sample1_summary, 
            "annotation/sample1_ChIP_annotation_summary.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)

if (exists("sample2_df")) {
    sample2_summary <- create_annotation_summary(sample2_df, "Sample2")
    write.table(sample2_summary, 
                "annotation/sample2_ChIP_annotation_summary.txt",
                sep = "\t", row.names = FALSE, quote = FALSE)
}

if (exists("sample3_df")) {
    sample3_summary <- create_annotation_summary(sample3_df, "Sample3")
    write.table(sample3_summary, 
                "annotation/sample3_ChIP_annotation_summary.txt",
                sep = "\t", row.names = FALSE, quote = FALSE)
}

cat("Annotation summaries created!\n")
```

## Step 6: Extract Target Genes

```r
# Extract genes associated with peaks for functional analysis
extract_target_genes <- function(annotation_df, sample_name) {
    
    # Filter for peaks associated with genes (exclude intergenic)
    gene_associated_peaks <- annotation_df[!grepl("Intergenic", annotation_df$annotation), ]
    
    # Extract unique gene IDs
    target_genes <- unique(gene_associated_peaks$geneId)
    target_genes <- target_genes[!is.na(target_genes) & target_genes != ""]
    
    # Save gene list
    write.table(target_genes, 
                paste0("functional_analysis/", sample_name, "_target_genes.txt"),
                row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    cat("Sample:", sample_name, "- Target genes:", length(target_genes), "\n")
    
    return(target_genes)
}

# Extract target genes for each sample
sample1_genes <- extract_target_genes(sample1_df, "sample1_ChIP")

if (exists("sample2_df")) {
    sample2_genes <- extract_target_genes(sample2_df, "sample2_ChIP")
}

if (exists("sample3_df")) {
    sample3_genes <- extract_target_genes(sample3_df, "sample3_ChIP")
}
```

## Step 7: Gene Ontology Enrichment Analysis

```r
# Function to perform GO enrichment analysis
perform_go_analysis <- function(gene_list, sample_name) {
    
    cat("Performing GO analysis for", length(gene_list), "genes from", sample_name, "\n")
    
    # Skip analysis if too few genes
    if (length(gene_list) < 10) {
        cat("Too few genes for meaningful analysis\n")
        return(NULL)
    }
    
    # Biological Process enrichment
    go_bp <- enrichGO(gene = gene_list,
                      OrgDb = org.Zm.eg.db,
                      ont = "BP",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.2,
                      readable = TRUE)
    
    # Save results if significant terms found
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
        
        # Save detailed results
        write.table(go_bp@result, 
                   paste0("functional_analysis/", sample_name, "_GO_BP.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        # Create visualization plots
        pdf(paste0("functional_analysis/", sample_name, "_GO_plots.pdf"), 
            width = 12, height = 8)
        
        # Bar plot of top terms
        print(barplot(go_bp, showCategory = 15, 
                     title = paste("Biological Processes -", sample_name)) +
              theme(axis.text.x = element_text(angle = 45, hjust = 1)))
        
        # Dot plot showing significance and gene ratio
        print(dotplot(go_bp, showCategory = 20,
                     title = paste("GO Enrichment -", sample_name)))
        
        dev.off()
        
        cat("GO analysis completed for", sample_name, "\n")
        return(go_bp)
    } else {
        cat("No significant GO terms found for", sample_name, "\n")
        return(NULL)
    }
}

# Perform GO analysis for each sample
sample1_go <- perform_go_analysis(sample1_genes, "sample1_ChIP")

if (exists("sample2_genes")) {
    sample2_go <- perform_go_analysis(sample2_genes, "sample2_ChIP")
}

if (exists("sample3_genes")) {
    sample3_go <- perform_go_analysis(sample3_genes, "sample3_ChIP")
}
```

## Step 8: Examine Specific Gene Categories

```r
# Look for auxin-related genes (since ARF27 is an auxin response factor)
find_auxin_genes <- function(annotation_df, sample_name) {
    
    # Look for genes with "auxin" in their description
    auxin_genes <- annotation_df[grepl("auxin|Auxin|AUXIN", 
                                       annotation_df$SYMBOL, 
                                       ignore.case = TRUE), ]
    
    if (nrow(auxin_genes) > 0) {
        cat("Found", nrow(auxin_genes), "auxin-related genes in", sample_name, "\n")
        
        # Save auxin-related genes
        write.table(auxin_genes[, c("seqnames", "start", "end", "geneId", "SYMBOL")], 
                   paste0("functional_analysis/", sample_name, "_auxin_genes.txt"),
                   sep = "\t", row.names = FALSE, quote = FALSE)
        
        return(auxin_genes)
    } else {
        cat("No obvious auxin-related genes found in", sample_name, "\n")
        return(NULL)
    }
}

# Search for auxin-related genes
sample1_auxin <- find_auxin_genes(sample1_df, "sample1_ChIP")

if (exists("sample2_df")) {
    sample2_auxin <- find_auxin_genes(sample2_df, "sample2_ChIP")
}

if (exists("sample3_df")) {
    sample3_auxin <- find_auxin_genes(sample3_df, "sample3_ChIP")
}
```

## Step 9: Peak Location Analysis

```r
# Analyze peak distribution across chromosomes
analyze_chromosome_distribution <- function(annotation_df, sample_name) {
    
    # Count peaks per chromosome
    chr_counts <- table(annotation_df$seqnames)
    chr_df <- data.frame(
        Chromosome = names(chr_counts),
        Peak_Count = as.numeric(chr_counts)
    )
    
    # Sort by chromosome number (handle both numeric and character names)
    chr_df$Chr_Num <- as.numeric(gsub("chr", "", chr_df$Chromosome))
    chr_df <- chr_df[order(chr_df$Chr_Num, na.last = TRUE), ]
    
    # Save chromosome distribution
    write.table(chr_df, 
               paste0("annotation/", sample_name, "_chromosome_distribution.txt"),
               sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("Chromosome distribution analysis completed for", sample_name, "\n")
    return(chr_df)
}

# Analyze chromosome distribution
sample1_chr <- analyze_chromosome_distribution(sample1_df, "sample1_ChIP")

if (exists("sample2_df")) {
    sample2_chr <- analyze_chromosome_distribution(sample2_df, "sample2_ChIP")
}

if (exists("sample3_df")) {
    sample3_chr <- analyze_chromosome_distribution(sample3_df, "sample3_ChIP")
}
```

## Step 10: Create Final Summary Report

```r
# Create comprehensive summary
create_final_summary <- function() {
    
    cat("Creating final annotation summary...\n")
    
    # Collect all summary information
    summary_text <- c(
        "ARF27 ChIP-seq ANNOTATION SUMMARY",
        "==================================",
        "",
        paste("Analysis Date:", Sys.Date()),
        "",
        "PEAK ANNOTATION RESULTS:",
        "------------------------"
    )
    
    # Add sample-specific summaries
    if (exists("sample1_df")) {
        promoter_peaks <- sum(grepl("Promoter", sample1_df$annotation))
        gene_body_peaks <- sum(grepl("Exon|Intron", sample1_df$annotation))
        intergenic_peaks <- sum(grepl("Intergenic", sample1_df$annotation))
        
        summary_text <- c(summary_text,
                         paste("Sample1 - Total peaks:", nrow(sample1_df)),
                         paste("  Promoter peaks:", promoter_peaks),
                         paste("  Gene body peaks:", gene_body_peaks),
                         paste("  Intergenic peaks:", intergenic_peaks),
                         paste("  Target genes:", length(sample1_genes)),
                         "")
    }
    
    # Add GO analysis summary
    summary_text <- c(summary_text,
                     "FUNCTIONAL ANALYSIS:",
                     "--------------------")
    
    if (exists("sample1_go") && !is.null(sample1_go)) {
        sig_terms <- nrow(sample1_go@result)
        summary_text <- c(summary_text,
                         paste("Sample1 - Significant GO terms:", sig_terms))
    }
    
    # Write summary to file
    writeLines(summary_text, "annotation/final_annotation_summary.txt")
    
    cat("Final summary created!\n")
}

# Create the final summary
create_final_summary()
```

## Step 11: Exit R and Return to Command Line

```r
# Save workspace for future reference
save.image("annotation/annotation_workspace.RData")

# Exit R
quit(save = "no")
```

## Step 12: Examine Results from Command Line

Now back in the terminal, let's examine what we created:

```bash
# Look at the annotation output files
ls annotation/
ls functional_analysis/

# Check the final summary
cat annotation/final_annotation_summary.txt

# Look at the first few annotated peaks for sample1
head -10 annotation/sample1_ChIP_detailed_annotation.txt

# Count different types of annotations
echo "Promoter peaks:"
grep -c "Promoter" annotation/sample1_ChIP_detailed_annotation.txt

echo "Gene body peaks:"
grep -c -E "(Exon|Intron)" annotation/sample1_ChIP_detailed_annotation.txt

echo "Intergenic peaks:"
grep -c "Intergenic" annotation/sample1_ChIP_detailed_annotation.txt
```

## Understanding Annotation Results

### Key Output Files

**detailed_annotation.txt**: Complete annotation information for each peak
- Genomic coordinates and nearby genes
- Distance to transcription start sites
- Gene symbols and descriptions
- Annotation category (promoter, exon, intron, intergenic)

**annotation_summary.txt**: Summary statistics
- Number of peaks in each category
- Percentage distribution of peak locations
- Total number of target genes

**GO_BP.txt**: Gene Ontology biological process enrichment
- Statistically overrepresented biological processes
- P-values and false discovery rates
- Gene lists for each significant term

**target_genes.txt**: Lists of genes associated with peaks
- Input for additional analyses
- Can be used for pathway analysis
- Represents direct regulatory targets

### Interpreting Peak Distribution

**High Promoter Enrichment** (>40% of peaks):
- Suggests direct transcriptional regulation
- Typical pattern for sequence-specific transcription factors
- Good for identifying primary regulatory targets

**Moderate Promoter Enrichment** (20-40%):
- Mix of direct and indirect regulation
- May include enhancer binding
- Still indicates significant regulatory activity

**Low Promoter Enrichment** (<20%):
- May indicate primarily enhancer binding
- Could suggest cooperative binding with other factors
- Might represent specific biological context

### Expected Results for ARF27

As an auxin response factor, you might expect:
- Enrichment in promoters of auxin-responsive genes
- GO terms related to:
  - Auxin signaling pathways
  - Root development
  - Gravitropism
  - Cell expansion
  - Hormone response

## Troubleshooting Common Issues

### No Significant GO Terms

```bash
# Check if you have enough target genes
wc -l functional_analysis/*_target_genes.txt

# Look at a few target genes to verify they're real
head -10 functional_analysis/sample1_ChIP_target_genes.txt
```

### R Package Installation Issues

```bash
# If R packages fail to install, try updating Bioconductor
R -e "BiocManager::install(version = 'devel')"
R -e "BiocManager::install('ChIPseeker', force = TRUE)"
```

### GTF File Problems

```bash
# Verify GTF file format
head -5 reference_genome/maize_annotations.gtf
grep -c "gene" reference_genome/maize_annotations.gtf
```

## Next Steps

With annotation complete, you now understand:
- Where ARF27 binds relative to genes
- Which genes are likely regulated by ARF27
- What biological processes ARF27 controls
- How ARF27 fits into cellular regulatory networks

This biological context transforms your peak coordinates into meaningful insights about ARF27's role in auxin signaling and plant development.

### Files for Further Analysis

- **Peak annotations**: For detailed investigation of specific targets
- **Target gene lists**: For additional pathway analysis or experimental validation
- **GO enrichment results**: For understanding biological significance
- **Visualization plots**: For presentations and publications

Your ARF27 ChIP-seq analysis is now complete, providing a comprehensive map of where this important transcription factor binds and what processes it regulates in maize.
