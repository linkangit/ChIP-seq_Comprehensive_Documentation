# Document 4: Reference Genome Preparation and Indexing

## Understanding Reference Genomes in ChIP-seq

The reference genome is like a detailed map that helps us determine where each sequencing read originated. For your ARF27 ChIP-seq analysis, we need both the DNA sequence of the maize genome and information about where genes are located within that sequence.

Think of the reference genome as a coordinate system. When we find that ARF27 binds to a specific location, we can use the reference genome to determine:
- Which chromosome the binding site is on
- The exact position within that chromosome
- Which genes are nearby
- Whether the binding site is in a promoter, gene body, or intergenic region

### Why Maize Genome Preparation is Challenging

The maize genome presents unique challenges compared to model organisms like Arabidopsis:

**Genome Size**: At approximately 2.3 billion base pairs, the maize genome is about 7 times larger than Arabidopsis and comparable to the human genome. This means:
- Longer processing times for all steps
- Higher memory requirements (16-32GB RAM recommended)
- Larger file sizes throughout the analysis
- More complex repetitive regions

**Repetitive Content**: About 85% of the maize genome consists of repetitive elements, primarily transposons. This creates challenges because:
- Many sequencing reads could potentially map to multiple locations
- Alignment algorithms must be carefully tuned to handle repetitive sequences
- Peak calling becomes more complex due to repetitive background

**Evolutionary History**: Maize is an ancient polyploid, meaning its genome contains duplicated regions from past genome doubling events. This results in:
- Gene families with multiple similar members
- Potential for reads to map to paralogous (similar but distinct) genes
- Complex interpretation of binding patterns

### The Two Essential Genome Components

For ChIP-seq analysis, you need two types of genomic information:

**1. Genome Sequence (FASTA file)**:
- Contains the actual DNA sequence of each chromosome
- Used by alignment programs to map reads to specific locations
- Serves as the reference coordinate system for all analyses

**2. Gene Annotations (GTF file)**:
- Describes where genes are located within the genome sequence
- Defines exon/intron boundaries, transcription start sites, and gene names
- Essential for determining which genes are near ARF27 binding sites
- Required for functional analysis and biological interpretation

## Setting Up the Reference Directory

Let's create a organized structure for our reference files:

```bash
# Make sure we're in the right place
cd chipseq_analysis
conda activate chipseq

# Create directory for reference genome files
mkdir -p reference_genome
cd reference_genome
```

## Step 1: Downloading the Maize Reference Genome

We'll use the latest high-quality maize reference genome (B73 RefGen v5) from Ensembl Plants, which provides well-annotated, standardized genomic data.

### Download the Genome Sequence

```bash
# Download the maize genome sequence (this will take 10-20 minutes)
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
```

**What this file contains**:
- All 10 maize chromosomes plus additional scaffolds
- Complete DNA sequence in FASTA format
- Compressed to save space (about 600MB compressed, 2.3GB uncompressed)

**Why we use this specific version**:
- B73 is the reference inbred line for maize genomics
- NAM-5.0 is the most recent, highest-quality assembly
- Ensembl provides consistent annotation across plant species

### Download the Gene Annotations

```bash
# Download the gene annotation file (this will take 2-5 minutes)
wget ftp://ftp.ensemblgenomes.org/pub/plants/release-54/gtf/zea_mays/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.54.gtf.gz
```

**What this file contains**:
- Locations of all known and predicted genes
- Exon and intron boundaries
- Transcription start and end sites
- Gene names, descriptions, and functional annotations
- Alternative transcript isoforms

### Decompress and Rename Files

```bash
# Decompress the downloaded files
gunzip Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
gunzip Zea_mays.Zm-B73-REFERENCE-NAM-5.0.54.gtf.gz

# Rename to simpler, more memorable names
mv Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa maize_genome.fa
mv Zea_mays.Zm-B73-REFERENCE-NAM-5.0.54.gtf maize_annotations.gtf
```

### Verify the Downloads

```bash
# Check that files were downloaded correctly
ls -lh *.fa *.gtf

# Look at the beginning of the genome file
head -5 maize_genome.fa

# Count how many chromosomes/scaffolds are in the genome
grep "^>" maize_genome.fa | wc -l

# Look at the annotation file structure
head -10 maize_annotations.gtf
```

**Expected output**:
- Genome file should be ~2.3GB
- Annotation file should be ~100-200MB
- You should see chromosome names like "1", "2", "3", etc.
- Annotation file should show gene coordinates and features

## Step 2: Understanding the File Formats

### FASTA Format Structure

The genome file uses FASTA format:
```
>1 dna:chromosome chromosome:Zm-B73-REFERENCE-NAM-5.0:1:1:308452471:1 REF
AACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC
CCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCT
```

**Header line** (starts with `>`): Contains chromosome name and metadata
**Sequence lines**: Actual DNA sequence using A, T, G, C, and N (unknown) bases

### GTF Format Structure

The annotation file uses GTF (Gene Transfer Format):
```
1	ensembl	gene	134	2764	.	+	.	gene_id "Zm00001eb000010"; gene_version "1"; gene_name "AC148152.3_FG001"; gene_source "ensembl"; gene_biotype "protein_coding";
1	ensembl	transcript	134	2764	.	+	.	gene_id "Zm00001eb000010"; transcript_id "Zm00001eb000010_T001"; 
```

**Columns**:
1. Chromosome name
2. Source (ensembl)
3. Feature type (gene, transcript, exon, etc.)
4. Start position
5. End position
6. Score (usually ".")
7. Strand (+ or -)
8. Frame (for coding sequences)
9. Attributes (gene names, IDs, etc.)

## Step 3: Creating the Bowtie2 Index

The genome index is a data structure that allows rapid searching of the genome sequence. Think of it like an index in a book - instead of reading every page to find a word, the index tells you exactly which pages contain that word.

### Why Indexing is Necessary

Without an index, aligning millions of reads would require:
1. Comparing each read to every position in the 2.3 billion base genome
2. This would take weeks or months of computation time
3. The process would be computationally impossible for large genomes

With an index:
1. The alignment algorithm can quickly identify candidate locations
2. Only promising locations need detailed comparison
3. Alignment completes in hours rather than weeks

### Building the Index

```bash
# Create directory for the index files
mkdir -p bowtie2_index

# Build the Bowtie2 index (this takes 1-3 hours and uses lots of memory)
bowtie2-build --threads 4 maize_genome.fa bowtie2_index/maize_genome
```

**Parameter explanation**:
- `--threads 4`: Use 4 CPU cores to speed up the process
- `maize_genome.fa`: Input genome file
- `bowtie2_index/maize_genome`: Prefix for output index files

**What happens during indexing**:
1. Bowtie2 reads the entire genome sequence
2. It creates multiple data structures for fast searching
3. Several index files are generated (6 files total)
4. The process requires 8-16GB of RAM for the maize genome

**Monitor progress**:
```bash
# Check how much memory is being used
free -h

# Monitor CPU usage
top

# Check if index files are being created
ls -lh bowtie2_index/
```

**Expected index files**:
```
maize_genome.1.bt2
maize_genome.2.bt2
maize_genome.3.bt2
maize_genome.4.bt2
maize_genome.rev.1.bt2
maize_genome.rev.2.bt2
```

### Verify Index Creation

```bash
# Check that all index files were created
ls -lh bowtie2_index/

# Verify index integrity (optional but recommended)
bowtie2-inspect -s bowtie2_index/maize_genome
```

## Step 4: Creating Additional Reference Files

### Genome Size File

Many tools need to know the size of each chromosome:

```bash
# Create FASTA index (generates .fai file)
samtools faidx maize_genome.fa

# Create chromosome sizes file
cut -f1,2 maize_genome.fa.fai > maize_genome.sizes

# Look at the chromosome sizes
head maize_genome.sizes
```

**What this creates**:
- `maize_genome.fa.fai`: Index file for rapid access to specific genome regions
- `maize_genome.sizes`: Two-column file with chromosome names and lengths

### Understanding Chromosome Naming

```bash
# See all chromosome/scaffold names
cut -f1 maize_genome.fa.fai | head -20

# Count chromosomes vs scaffolds
echo "Main chromosomes:"
cut -f1 maize_genome.fa.fai | grep -E "^[0-9]+$" | wc -l

echo "Additional scaffolds:"
cut -f1 maize_genome.fa.fai | grep -v -E "^[0-9]+$" | wc -l
```

**Expected output**:
- 10 main chromosomes (named 1, 2, 3, ..., 10)
- Additional scaffolds and organellar genomes (chloroplast, mitochondria)

## Step 5: Validating Your Reference Setup

### Check File Integrity

```bash
# Verify genome file is complete
tail -5 maize_genome.fa

# Check annotation file format
grep -c "^#" maize_annotations.gtf  # Should be 0 (no header lines)
grep -c "gene" maize_annotations.gtf  # Should be ~40,000-50,000

# Verify index files exist and have reasonable sizes
ls -lh bowtie2_index/*.bt2
```

### Test the Index

```bash
# Quick test alignment (should complete without errors)
echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCG" > test_read.fa
bowtie2 -x bowtie2_index/maize_genome -f test_read.fa -S test_alignment.sam

# Clean up test files
rm test_read.fa test_alignment.sam
```

### Memory and Disk Space Check

```bash
# Check disk space usage
du -sh bowtie2_index/
du -sh maize_genome.*
du -sh maize_annotations.gtf

# Verify you have enough disk space for analysis
df -h .
```

**Expected disk usage**:
- Genome file: ~2.3GB
- Annotation file: ~200MB
- Bowtie2 index: ~4-6GB
- Total: ~7-9GB for reference files

## Understanding Index Parameters and Optimization

### Why We Use These Bowtie2 Settings

**Default indexing parameters work well for maize because**:
- Bowtie2 automatically optimizes for genome size
- Default settings balance memory usage vs. speed
- The resulting index handles repetitive sequences appropriately

**When you might need different settings**:
- Very limited RAM: Use `bowtie2-build --bmax 40000000` to reduce memory
- Faster indexing: Add more threads with `--threads 8` if you have more cores
- Different genome: These same steps work for any plant genome

### Index Performance Characteristics

The Bowtie2 index provides:
- **Fast alignment**: Can align 1-2 million reads per minute
- **Memory efficiency**: Index fits in 4-6GB RAM during alignment
- **Accuracy**: Handles up to 2-3 mismatches per read
- **Repeat handling**: Manages repetitive sequences appropriately for ChIP-seq

## Troubleshooting Common Issues

### Download Problems

```bash
# If download fails, check internet connection and retry
wget --continue ftp://ftp.ensemblgenomes.org/pub/plants/release-54/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# If server is slow, try alternative mirror
wget http://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-54/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz
```

### Memory Issues During Indexing

```bash
# If indexing fails due to memory, try reducing threads
bowtie2-build --threads 2 maize_genome.fa bowtie2_index/maize_genome

# Or reduce memory usage
bowtie2-build --bmax 40000000 --threads 4 maize_genome.fa bowtie2_index/maize_genome
```

### Disk Space Issues

```bash
# Check available space
df -h

# Clean up unnecessary files if needed
rm *.gz  # Remove compressed files after extraction

# Use symbolic links if files are stored elsewhere
ln -s /path/to/large/storage/maize_genome.fa .
```

## Preparing for Alignment

After completing this step, you should have:

### Required Files Checklist

```bash
# Verify all required files exist
ls -la maize_genome.fa          # Genome sequence
ls -la maize_annotations.gtf    # Gene annotations  
ls -la maize_genome.fa.fai      # Genome index
ls -la maize_genome.sizes       # Chromosome sizes
ls -la bowtie2_index/*.bt2      # Bowtie2 index files (6 files)
```

### Quality Checks

```bash
# File size verification
echo "Genome file size (should be ~2.3GB):"
ls -lh maize_genome.fa

echo "Number of chromosomes (should be ~10 + scaffolds):"
grep -c "^>" maize_genome.fa

echo "Number of genes (should be ~40,000-50,000):"
grep -c "gene" maize_annotations.gtf

echo "Index files (should be 6 files, 4-6GB total):"
ls -lh bowtie2_index/
```

### Documentation

```bash
# Create a reference setup summary
cat > reference_setup_summary.txt << 'EOF'
REFERENCE GENOME SETUP SUMMARY
==============================

Date: 
Genome version: Zea mays B73 RefGen v5 (NAM-5.0)
Annotation version: Ensembl Plants release 54
Source: ftp://ftp.ensemblgenomes.org/

Files created:
- maize_genome.fa (genome sequence)
- maize_annotations.gtf (gene annotations)
- maize_genome.fa.fai (sequence index)
- maize_genome.sizes (chromosome sizes)
- bowtie2_index/ (alignment index)

Index parameters:
- bowtie2-build with 4 threads
- Default memory settings
- Standard parameters for large plant genome

Verification:
- All files present and correct sizes
- Index integrity verified
- Test alignment successful
EOF

# Add file sizes to summary
echo "" >> reference_setup_summary.txt
echo "File sizes:" >> reference_setup_summary.txt
ls -lh maize_genome.fa maize_annotations.gtf >> reference_setup_summary.txt
echo "Index size:" >> reference_setup_summary.txt
du -sh bowtie2_index/ >> reference_setup_summary.txt
```

## Next Steps

With your reference genome properly prepared and indexed, you're now ready to align your processed ChIP-seq reads to identify where they originated in the maize genome. The index you've created will make this alignment process efficient and accurate.

The reference setup is a one-time process - you can reuse these same reference files for future maize ChIP-seq experiments, saving significant time and computational resources.

### Return to Analysis Directory

```bash
# Go back to the main analysis directory for the next step
cd ..
pwd  # Should show chipseq_analysis
```

Your reference genome is now ready for read alignment, which will reveal the genomic locations where ARF27 binds in maize.
