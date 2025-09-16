# Step 2: Read Mapping to Maize B73-v4 Reference Genome

## Objective
Align quality-controlled ChIP-seq reads to the maize B73-v4 reference genome to determine the genomic locations where ARF27 binds.

## Why Read Mapping is Critical

Accurate read mapping is essential for ChIP-seq analysis because:
- It determines the precise genomic coordinates of potential binding sites
- Poor mapping leads to false positive and false negative peaks
- Mapping quality affects downstream peak calling sensitivity
- ARF27 binding sites need to be mapped accurately to identify target genes

## Understanding the Maize B73-v4 Reference Genome

The B73-v4 reference genome features:
- **Size**: ~2.1 Gb with 10 chromosomes
- **Gene models**: ~39,000 protein-coding genes
- **Repetitive content**: ~85% (high compared to other plant genomes)
- **Chromosome naming**: chr1, chr2, ..., chr10
- **Coordinate system**: 1-based positioning

## Mapping Strategy for ARF27 ChIP-seq

### Why BWA is Recommended:
- Excellent performance with short reads (50-150 bp)
- Good handling of repetitive sequences (important for maize)
- Widely used and well-validated for ChIP-seq
- Fast and memory-efficient

### Alternative: Bowtie2
- Also excellent for ChIP-seq
- Sometimes better for longer reads
- Good sensitivity for gapped alignments

## Step-by-Step Instructions

### 1. Prepare the Reference Genome

First, download and prepare the maize B73-v4 reference genome:

```bash
# Navigate to a reference directory (create if needed)
mkdir -p ../reference_genome
cd ../reference_genome

# Download maize B73-v4 reference genome
# Option 1: From MaizeGDB
wget https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Option 2: From Ensembl Plants
# wget http://ftp.ensemblgenomes.org/pub/plants/release-52/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz

# Extract the genome file
gunzip Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Rename for clarity
mv Zm-B73-REFERENCE-NAM-5.0.fa maize_B73_v4.fa

# Check genome file
head -20 maize_B73_v4.fa
grep ">" maize_B73_v4.fa | head -15
```

### 2. Install and Set Up BWA

```bash
# Install BWA using conda (recommended)
conda install -c bioconda bwa

# Or install manually from: http://bio-bwa.sourceforge.net/
```

### 3. Index the Reference Genome

This step creates index files that BWA uses for fast alignment:

```bash
# Index the maize reference genome
# This takes 30-60 minutes and ~6GB RAM
bwa index maize_B73_v4.fa

# Check that index files were created
ls -la maize_B73_v4.fa*
```

**Expected output files:**
- `maize_B73_v4.fa.amb`
- `maize_B73_v4.fa.ann`
- `maize_B73_v4.fa.bwt`
- `maize_B73_v4.fa.pac`
- `maize_B73_v4.fa.sa`

### 4. Prepare for Mapping

Navigate back to your main analysis directory:

```bash
cd ../chipseq_arf27_analysis
```

### 5. Map Reads with Quality Trimming (if needed)

Based on your QC results from Step 1, choose the appropriate mapping approach:

#### Option A: Direct Mapping (if QC was good)

```bash
# Map ARF27 ChIP replicate 1
echo "Mapping ARF27 ChIP replicate 1..."
bwa mem -t 8 \
    -M \
    ../reference_genome/maize_B73_v4.fa \
    01_raw_data/ARF27_ChIP_rep1_R1.fastq.gz \
    01_raw_data/ARF27_ChIP_rep1_R2.fastq.gz \
    > 04_mapped_data/ARF27_ChIP_rep1.sam 2> logs/ARF27_ChIP_rep1_mapping.log

# Map ARF27 ChIP replicate 2
echo "Mapping ARF27 ChIP replicate 2..."
bwa mem -t 8 \
    -M \
    ../reference_genome/maize_B73_v4.fa \
    01_raw_data/ARF27_ChIP_rep2_R1.fastq.gz \
    01_raw_data/ARF27_ChIP_rep2_R2.fastq.gz \
    > 04_mapped_data/ARF27_ChIP_rep2.sam 2> logs/ARF27_ChIP_rep2_mapping.log

# Map Input control replicate 1
echo "Mapping Input control replicate 1..."
bwa mem -t 8 \
    -M \
    ../reference_genome/maize_B73_v4.fa \
    01_raw_data/Input_control_rep1_R1.fastq.gz \
    01_raw_data/Input_control_rep1_R2.fastq.gz \
    > 04_mapped_data/Input_control_rep1.sam 2> logs/Input_control_rep1_mapping.log

# Map Input control replicate 2
echo "Mapping Input control replicate 2..."
bwa mem -t 8 \
    -M \
    ../reference_genome/maize_B73_v4.fa \
    01_raw_data/Input_control_rep2_R1.fastq.gz \
    01_raw_data/Input_control_rep2_R2.fastq.gz \
    > 04_mapped_data/Input_control_rep2.sam 2> logs/Input_control_rep2_mapping.log
```

#### Option B: Mapping with Adapter Trimming (if adapters detected)

If FastQC revealed adapter contamination, use this approach:

```bash
# Install Trimmomatic if not available
conda install -c bioconda trimmomatic

# Create adapter file (adjust sequences based on your adapters)
cat > adapters.fa << EOF
>TruSeq3-PE-2
AGATCGGAAGAGCGGTTCAGCAGGAATGCCGAGACCGATCTCGTATGCCGTCTTCTGCTTG
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
EOF

# Trim and map ARF27 ChIP replicate 1
echo "Trimming and mapping ARF27 ChIP replicate 1..."
trimmomatic PE -phred33 \
    01_raw_data/ARF27_ChIP_rep1_R1.fastq.gz \
    01_raw_data/ARF27_ChIP_rep1_R2.fastq.gz \
    03_trimmed_data/ARF27_ChIP_rep1_R1_paired.fastq.gz \
    03_trimmed_data/ARF27_ChIP_rep1_R1_unpaired.fastq.gz \
    03_trimmed_data/ARF27_ChIP_rep1_R2_paired.fastq.gz \
    03_trimmed_data/ARF27_ChIP_rep1_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30

# Map trimmed reads
bwa mem -t 8 \
    -M \
    ../reference_genome/maize_B73_v4.fa \
    03_trimmed_data/ARF27_ChIP_rep1_R1_paired.fastq.gz \
    03_trimmed_data/ARF27_ChIP_rep1_R2_paired.fastq.gz \
    > 04_mapped_data/ARF27_ChIP_rep1.sam 2> logs/ARF27_ChIP_rep1_mapping.log

# Repeat for all other samples...
```

**BWA Command Explanation:**
- `-t 8`: Use 8 CPU threads (adjust based on your system)
- `-M`: Mark shorter split hits as secondary (recommended for Picard compatibility)
- Input order: reference genome, forward reads, reverse reads
- `2>`: Redirect mapping statistics to log file

**Trimmomatic Parameters:**
- `ILLUMINACLIP:adapters.fa:2:30:10`: Remove adapters
- `LEADING:3`: Cut bases with quality <3 from start
- `TRAILING:3`: Cut bases with quality <3 from end
- `SLIDINGWINDOW:4:15`: Cut when 4-base window drops below quality 15
- `MINLEN:30`: Drop reads shorter than 30 bp

### 6. Check Mapping Progress and Statistics

Monitor mapping progress:

```bash
# Check if mapping is complete
ls -la 04_mapped_data/

# Quick check of mapping log
tail -20 logs/ARF27_ChIP_rep1_mapping.log

# Count mapped reads in SAM file
grep -v "^@" 04_mapped_data/ARF27_ChIP_rep1.sam | wc -l
```

## Interpreting Mapping Statistics

Examine the mapping log files for these key metrics:

### Good Mapping Indicators:
- **Overall alignment rate**: >85% for high-quality data
- **Properly paired**: >80% for paired-end data
- **Mapping quality**: Most reads with MAPQ ≥ 10

### Warning Signs:
- **Low alignment rate**: <70% (may indicate quality issues or wrong reference)
- **High mismatch rate**: >5% (check for contamination)
- **Very low properly paired**: <60% (fragment size issues)

### Expected Results for Maize ChIP-seq:
```
Mapping Statistics Example:
- Total reads: 45,234,892
- Mapped reads: 39,450,234 (87.2%)
- Properly paired: 37,123,456 (82.1%)
- Singleton reads: 2,326,778 (5.1%)
```

## Quality Control of Mapped Data

### 1. Basic Mapping Statistics

```bash
# Install samtools if not available
conda install -c bioconda samtools

# Get basic statistics for each sample
for sample in ARF27_ChIP_rep1 ARF27_ChIP_rep2 Input_control_rep1 Input_control_rep2
do
    echo "Statistics for $sample:"
    samtools flagstat 04_mapped_data/${sample}.sam
    echo "----------------------------------------"
done
```

### 2. Check Insert Size Distribution

For paired-end ChIP-seq, insert sizes should typically be 150-300 bp:

```bash
# Convert SAM to BAM for insert size analysis
samtools view -bS 04_mapped_data/ARF27_ChIP_rep1.sam > temp_rep1.bam
samtools sort temp_rep1.bam -o temp_rep1_sorted.bam

# Generate insert size statistics
samtools stats temp_rep1_sorted.bam | grep "insert size"

# Clean up temporary files
rm temp_rep1.bam temp_rep1_sorted.bam
```

## Common Issues and Troubleshooting

### Issue 1: Low Mapping Rate (<70%)

**Possible Causes:**
- Wrong reference genome
- Poor read quality
- Significant adapter contamination
- Sample contamination

**Diagnosis Commands:**
```bash
# Check unmapped reads
samtools view -f 4 04_mapped_data/ARF27_ChIP_rep1.sam | head -10

# BLAST a few unmapped reads to identify source
# Extract unmapped read sequences and BLAST against NCBI
```

**Solutions:**
- Verify you're using the correct B73-v4 reference
- Increase trimming stringency
- Check for contamination

### Issue 2: Very High Mapping Rate (>98%)

**Possible Concerns:**
- May indicate over-trimming
- Verify this is consistent across samples
- Check mapping quality distribution

### Issue 3: SAM File Size Issues

**Large SAM files (>50GB) per sample:**
- Normal for high-coverage ChIP-seq
- Ensure sufficient disk space
- Consider streaming to BAM conversion

### Issue 4: Mapping Takes Too Long

**Optimization:**
```bash
# Increase threads (if you have more CPUs)
bwa mem -t 16 ...  # instead of -t 8

# Check system resources
htop
df -h  # Check disk space
```

## File Size Expectations

For typical ChIP-seq data:
- **Input FASTQ**: 2-8 GB per file
- **SAM files**: 20-80 GB per sample
- **Mapping time**: 2-6 hours per sample (depending on coverage and system)

## Analysis Log Update

Document your mapping results:

```bash
# Update analysis log
cat >> logs/step2_mapping_log.txt << EOF

Step 2: Read Mapping Results
Date: $(date)

Mapping Summary:
- Reference genome: maize B73-v4
- Aligner: BWA mem
- Trimming applied: [Yes/No - specify parameters]

Mapping Statistics:
- ARF27_ChIP_rep1: [%mapped] - [total reads]
- ARF27_ChIP_rep2: [%mapped] - [total reads]
- Input_control_rep1: [%mapped] - [total reads]
- Input_control_rep2: [%mapped] - [total reads]

Issues Encountered:
- [List any problems and solutions]

Quality Assessment:
- Overall mapping quality: [Good/Acceptable/Poor]
- Proceed to BAM conversion: [Yes/No]

EOF
```

## File Cleanup and Organization

```bash
# Check final file structure
ls -la 04_mapped_data/
ls -la logs/

# Optional: Compress SAM files if space is limited
# (Only after confirming quality)
# gzip 04_mapped_data/*.sam
```

## Next Steps Preparation

Before proceeding to Step 3 (SAM to BAM conversion), ensure:

- [ ] All samples mapped successfully
- [ ] Mapping rates are acceptable (>70%)
- [ ] Log files contain complete statistics
- [ ] Sufficient disk space for BAM files
- [ ] Documented any issues and solutions

## Key Takeaways

- Accurate mapping is foundation for all downstream analysis
- Maize's repetitive genome requires careful attention to mapping quality
- Document all parameters and statistics for reproducibility
- Quality control at each step prevents problems downstream

Ready for **Step 3: SAM to BAM Conversion**? We'll convert your aligned reads to the more efficient BAM format and prepare them for duplicate removal.
