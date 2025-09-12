# Document 5: Read Alignment to the Genome

## Understanding Read Alignment in ChIP-seq

Read alignment is the process of determining where each sequencing read originated in the genome. Think of it as solving millions of puzzles - for each short DNA sequence (read), we need to find its correct position among the 2.3 billion base pairs of the maize genome.

This step is crucial because it transforms abstract sequencing reads into genomic coordinates, allowing us to identify specific locations where ARF27 binds to DNA.

### The Biological Basis of ChIP-seq Alignment

In your ARF27 ChIP-seq experiment, each sequencing read represents:
- A fragment of DNA that was bound by ARF27 in living plant cells
- One end of a larger DNA fragment (typically 100-300 base pairs)
- A piece of evidence for transcription factor binding at that genomic location

When we align these reads to the genome, we're essentially reconstructing the map of where ARF27 was bound when the cells were cross-linked.

### Why Alignment is Challenging for Maize

**Genome Size**: With 2.3 billion base pairs, there are many potential locations where a short read might map. A 50-base read could theoretically match millions of locations by chance alone.

**Repetitive Sequences**: About 85% of the maize genome consists of repetitive elements. This means:
- Many reads could map to multiple locations with equal confidence
- We need strategies to handle "multi-mapping" reads appropriately
- Alignment parameters must balance sensitivity with specificity

**Sequencing Errors**: Even high-quality reads contain 1-2 errors per 100 bases. The alignment algorithm must:
- Allow for mismatches due to sequencing errors
- Distinguish between real genetic variation and technical errors
- Maintain accuracy while being tolerant of minor differences

## Bowtie2: The Alignment Algorithm

Bowtie2 is specifically designed for aligning short sequencing reads to large genomes. It uses sophisticated algorithms to quickly identify potential alignment locations and then carefully evaluate the best match for each read.

### How Bowtie2 Works

1. **Seed Finding**: Bowtie2 extracts short "seeds" from each read and uses the genome index to find potential alignment locations
2. **Extension**: For each promising location, it attempts to align the entire read
3. **Scoring**: It calculates alignment scores based on matches, mismatches, and gaps
4. **Selection**: It chooses the best alignment location for each read

### ChIP-seq Specific Alignment Considerations

**Fragment Size Awareness**: ChIP-seq reads come from larger DNA fragments. Bowtie2 can use this information to improve alignment accuracy.

**Quality-Based Alignment**: Higher quality bases are weighted more heavily in alignment decisions.

**Repetitive Region Handling**: Bowtie2 can report multiple alignments or randomly assign reads to one of several equally good locations.

## Setting Up for Alignment

```bash
# Make sure we're in the right directory and environment
cd chipseq_analysis
conda activate chipseq

# Create directory for alignment results
mkdir -p alignment
cd alignment
```

## Step 1: Understanding Your Processed Reads

Before alignment, let's examine what files we have from the trimming step:

```bash
# Look at what trimmed files are available
ls ../processed_data/

# Count reads in a processed file (example with first file)
first_file=$(ls ../processed_data/*_trimmed.fastq.gz | head -1)
echo "Counting reads in: $first_file"
total_lines=$(zcat $first_file | wc -l)
read_count=$((total_lines / 4))
echo "Total reads: $read_count"
```

## Step 2: Single-End Read Alignment

For single-end ChIP-seq data (most common), align each sample individually:

### Align First Sample

```bash
# Example: Align the first sample (replace sample1 with your actual sample name)
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -U ../processed_data/sample1_trimmed.fastq.gz \
    --threads 4 \
    --very-sensitive \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample1_sorted.bam -

# Index the BAM file for fast access
samtools index sample1_sorted.bam
```

### Align Additional Samples

```bash
# Align sample 2
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -U ../processed_data/sample2_trimmed.fastq.gz \
    --threads 4 \
    --very-sensitive \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample2_sorted.bam -

samtools index sample2_sorted.bam

# Align sample 3
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -U ../processed_data/sample3_trimmed.fastq.gz \
    --threads 4 \
    --very-sensitive \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample3_sorted.bam -

samtools index sample3_sorted.bam

# Continue for all your samples...
```

## Step 3: Paired-End Read Alignment (if applicable)

If your data is paired-end, align both read files together:

```bash
# Example: Align paired-end sample 1
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -1 ../processed_data/sample1_R1_paired.fastq.gz \
    -2 ../processed_data/sample1_R2_paired.fastq.gz \
    --threads 4 \
    --very-sensitive \
    -X 2000 \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample1_paired_sorted.bam -

samtools index sample1_paired_sorted.bam

# Align paired-end sample 2
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -1 ../processed_data/sample2_R1_paired.fastq.gz \
    -2 ../processed_data/sample2_R2_paired.fastq.gz \
    --threads 4 \
    --very-sensitive \
    -X 2000 \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample2_paired_sorted.bam -

samtools index sample2_paired_sorted.bam
```

## Understanding Alignment Parameters

Let's break down each parameter and understand why it's chosen for ChIP-seq:

### Core Alignment Parameters

**-x ../reference_genome/bowtie2_index/maize_genome**
- Specifies the path to the genome index we built in Document 4
- This tells Bowtie2 where to find the pre-computed index files

**-U ../processed_data/sample1_trimmed.fastq.gz**
- Input file for single-end reads
- Uses the trimmed, high-quality reads from Document 3

**-1 and -2** (for paired-end)
- -1: First read of paired-end data (R1 file)
- -2: Second read of paired-end data (R2 file)
- Bowtie2 aligns both reads together, improving accuracy

### Performance Parameters

**--threads 4**
- Uses 4 CPU cores for parallel processing
- Speeds up alignment significantly (can use 6-8 threads if available)
- More threads = faster alignment, but diminishing returns beyond 8

**-q**
- Indicates input files are in FASTQ format (with quality scores)
- Allows Bowtie2 to use quality information for better alignment decisions

### Sensitivity Parameters

**--very-sensitive**
- Uses thorough alignment parameters optimized for accuracy
- Equivalent to: `-D 20 -R 3 -N 0 -L 20 -i S,1,0.50`
- Slower than default but more accurate for ChIP-seq
- Essential for transcription factor ChIP-seq where precision matters

**Alternative sensitivity levels**:
- `--fast`: Faster but less accurate
- `--sensitive`: Good balance (default)
- `--very-sensitive`: Most accurate (recommended for ChIP-seq)
- `--very-sensitive-local`: For reads with adapters or poor ends

### Paired-End Specific Parameters

**-X 2000**
- Maximum insert size for paired-end alignment
- Sets the expected distance between read pairs
- 2000 bp accommodates most ChIP-seq fragment sizes
- Prevents alignment of reads that are too far apart

### SAMtools Pipeline Parameters

**samtools view -bS -**
- Converts SAM format to compressed BAM format
- `-b`: Output in BAM format
- `-S`: Input is SAM format
- `-`: Read from standard input (piped from Bowtie2)

**samtools sort -o sample1_sorted.bam -**
- Sorts alignments by genomic coordinate
- Sorted BAM files are required for most downstream tools
- `-o`: Specify output file name
- `-`: Read from standard input

## Step 4: Monitor Alignment Progress

### Check Alignment Status

```bash
# See what BAM files have been created
ls -lh *.bam

# Monitor alignment progress (while running)
top

# Check memory usage
free -h
```

### Understanding Alignment Output

During alignment, you'll see output like this:

```
25000000 reads; of these:
  25000000 (100.00%) were unpaired; of these:
    2341234 (9.36%) aligned 0 times
    20167891 (80.67%) aligned exactly 1 time
    2490875 (9.97%) aligned >1 times
90.64% overall alignment rate
```

**Key metrics to evaluate**:
- **Overall alignment rate**: Should be >70% for good ChIP-seq data
- **Aligned exactly 1 time**: Uniquely mapped reads (ideal for peak calling)
- **Aligned >1 times**: Multi-mapping reads (common in repetitive genomes)
- **Aligned 0 times**: Unmapped reads (contamination, low quality, or novel sequences)

## Step 5: Generate Alignment Statistics

### Create Detailed Statistics for Each Sample

```bash
# Generate statistics for sample 1
samtools flagstat sample1_sorted.bam > sample1_alignment_stats.txt

# Generate statistics for sample 2  
samtools flagstat sample2_sorted.bam > sample2_alignment_stats.txt

# Generate statistics for sample 3
samtools flagstat sample3_sorted.bam > sample3_alignment_stats.txt

# Continue for all samples...
```

### View Alignment Statistics

```bash
# Look at alignment statistics for first sample
cat sample1_alignment_stats.txt

# Create summary of all samples
echo "Sample\tTotal_Reads\tMapped_Reads\tMapping_Rate" > alignment_summary.txt

for bam in *.bam; do
    sample=$(basename $bam _sorted.bam)
    
    # Extract statistics from flagstat output
    total=$(samtools flagstat $bam | head -1 | cut -d' ' -f1)
    mapped=$(samtools flagstat $bam | head -5 | tail -1 | cut -d' ' -f1)
    rate=$(echo "scale=2; $mapped * 100 / $total" | bc)
    
    echo -e "$sample\t$total\t$mapped\t$rate%" >> alignment_summary.txt
done

# View the summary
column -t alignment_summary.txt
```

## Step 6: Post-Alignment Processing

### Remove PCR Duplicates

PCR duplicates are identical reads that arise from amplification artifacts rather than independent sequencing events. They can artificially inflate peak signals.

```bash
# Remove duplicates from sample 1
samtools rmdup sample1_sorted.bam sample1_dedup.bam
samtools index sample1_dedup.bam

# Remove duplicates from sample 2
samtools rmdup sample2_sorted.bam sample2_dedup.bam
samtools index sample2_dedup.bam

# Remove duplicates from sample 3
samtools rmdup sample3_sorted.bam sample3_dedup.bam
samtools index sample3_dedup.bam

# Continue for all samples...
```

### Filter for High-Quality Alignments

Remove low-quality and ambiguous alignments:

```bash
# Filter sample 1 for high-quality alignments
samtools view -b -q 20 -F 1024 sample1_dedup.bam > sample1_filtered.bam
samtools index sample1_filtered.bam

# Filter sample 2
samtools view -b -q 20 -F 1024 sample2_dedup.bam > sample2_filtered.bam
samtools index sample2_filtered.bam

# Filter sample 3
samtools view -b -q 20 -F 1024 sample3_dedup.bam > sample3_filtered.bam
samtools index sample3_filtered.bam

# Continue for all samples...
```

**Filter parameters explained**:
- `-q 20`: Keep only reads with mapping quality ≥20 (99% confidence)
- `-F 1024`: Exclude PCR duplicates
- `-b`: Output in BAM format

### Generate Final Statistics

```bash
# Statistics after filtering for each sample
samtools flagstat sample1_filtered.bam > sample1_final_stats.txt
samtools flagstat sample2_filtered.bam > sample2_final_stats.txt
samtools flagstat sample3_filtered.bam > sample3_final_stats.txt

# Create final summary
echo "Sample\tOriginal_Reads\tAfter_Dedup\tAfter_Filter\tFinal_Rate" > final_alignment_summary.txt

for sample in sample1 sample2 sample3; do
    original=$(samtools flagstat ${sample}_sorted.bam | head -1 | cut -d' ' -f1)
    dedup=$(samtools flagstat ${sample}_dedup.bam | head -1 | cut -d' ' -f1)
    final=$(samtools flagstat ${sample}_filtered.bam | head -1 | cut -d' ' -f1)
    rate=$(echo "scale=2; $final * 100 / $original" | bc)
    
    echo -e "$sample\t$original\t$dedup\t$final\t$rate%" >> final_alignment_summary.txt
done

# View final summary
column -t final_alignment_summary.txt
```

## Understanding Alignment Quality Metrics

### Mapping Quality Scores

Mapping quality (MAPQ) represents confidence in alignment position:
- **MAPQ 0**: Read maps to multiple locations with equal confidence
- **MAPQ 1-19**: Low confidence alignment
- **MAPQ 20-39**: Good confidence (99-99.9% probability of correct mapping)
- **MAPQ 40+**: High confidence (>99.9% probability)

### Expected Results for Good ChIP-seq Data

**Mapping Rates**:
- **>90%**: Excellent quality
- **80-90%**: Good quality  
- **70-80%**: Acceptable quality
- **<70%**: Poor quality, investigate issues

**Duplicate Rates**:
- **20-40%**: Normal for ChIP-seq
- **40-60%**: Acceptable, moderate enrichment
- **>60%**: High, may indicate over-amplification

**Final Read Retention**:
- **>60%**: Excellent after all filtering
- **40-60%**: Good retention
- **<40%**: Concerning, may indicate quality issues

## Troubleshooting Common Alignment Issues

### Low Mapping Rates

**Possible causes**:
- Adapter contamination not fully removed
- Wrong reference genome
- Poor sequencing quality
- Species contamination

**Solutions**:
```bash
# Check for remaining adapters
samtools view sample1_sorted.bam | head -1000 | cut -f10 | grep -i AGATCGG

# Try less stringent alignment
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -U ../processed_data/sample1_trimmed.fastq.gz \
    --threads 4 \
    --sensitive \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample1_lenient.bam -
```

### High Duplicate Rates

```bash
# Check duplication levels
samtools flagstat sample1_dedup.bam
samtools flagstat sample1_sorted.bam

# Calculate duplication rate
original=$(samtools flagstat sample1_sorted.bam | head -1 | cut -d' ' -f1)
dedup=$(samtools flagstat sample1_dedup.bam | head -1 | cut -d' ' -f1)
dup_rate=$(echo "scale=2; ($original - $dedup) * 100 / $original" | bc)
echo "Duplication rate: $dup_rate%"
```

### Memory Issues

```bash
# If running out of memory, reduce threads
bowtie2 -x ../reference_genome/bowtie2_index/maize_genome \
    -U ../processed_data/sample1_trimmed.fastq.gz \
    --threads 2 \
    --very-sensitive \
    -q \
    | samtools view -bS - \
    | samtools sort -o sample1_sorted.bam -
```

## Validating Alignment Results

### Quick Visual Check

```bash
# Look at a few aligned reads
samtools view sample1_filtered.bam | head -5

# Check alignment distribution across chromosomes
samtools idxstats sample1_filtered.bam | head -15
```

### Insert Size Analysis (for paired-end data)

```bash
# Generate insert size statistics for paired-end data
samtools stats sample1_paired_sorted.bam | grep "^IS" | cut -f2- > sample1_insert_sizes.txt

# View insert size distribution
head -20 sample1_insert_sizes.txt
```

## Preparing for Peak Calling

### Create Clean Working Directory

```bash
# Return to main analysis directory
cd ..

# Verify we have all the filtered BAM files we need
ls alignment/*_filtered.bam

# Check that all BAM files are indexed
ls alignment/*_filtered.bam.bai
```

### Documentation

```bash
# Create alignment summary document
cat > alignment/alignment_summary.txt << 'EOF'
ALIGNMENT PROCESS SUMMARY
========================

Date: 
Aligner: Bowtie2
Reference: Maize B73 RefGen v5
Parameters: --very-sensitive, -q 20, removed duplicates

Alignment Statistics:
EOF

# Add statistics for each sample
for bam in alignment/*_filtered.bam; do
    sample=$(basename $bam _filtered.bam)
    stats=$(samtools flagstat $bam | head -1)
    echo "$sample: $stats" >> alignment/alignment_summary.txt
done
```

### Files Ready for Peak Calling

You should now have:
- **Filtered BAM files**: High-quality, deduplicated alignments
- **BAM index files**: For rapid access to alignment data
- **Alignment statistics**: Documentation of mapping success
- **Quality metrics**: Evidence of successful alignment

These filtered BAM files contain the precise genomic locations where your sequencing reads aligned, representing the locations where ARF27 was bound to DNA in your experimental samples.

The next step will be peak calling, where we identify genomic regions with statistically significant enrichment of reads, indicating likely ARF27 binding sites.
