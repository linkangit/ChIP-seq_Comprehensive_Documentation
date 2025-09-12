# Document 3: Read Processing and Quality Improvement

## The Science Behind Read Trimming

Read processing transforms raw, imperfect sequencing data into clean, high-quality reads suitable for accurate genome alignment. Think of this step as data refinement - you're removing technical artifacts while preserving the biological signal that reveals where ARF27 binds in the maize genome.

### Why Raw Reads Need Processing

Sequencing technologies, while powerful, introduce several types of technical artifacts:

**Sequencing Quality Decline**: DNA polymerases become less efficient over longer read lengths, leading to increased error rates toward the 3' end of reads. This is particularly pronounced in Illumina sequencing due to:
- Accumulation of fluorescent dye molecules that interfere with imaging
- Progressive damage to DNA templates during repeated chemistry cycles
- Optical crowding on the flowcell surface

**Adapter Contamination**: Sequencing adapters are artificial DNA sequences added during library preparation. When the actual DNA fragment is shorter than the sequencing read length, the sequencer continues reading into the adapter sequence. This creates several problems:
- Adapter sequences don't exist in the genome, causing alignment failures
- Partial adapter sequences may align to spurious genomic locations
- Adapter dimers (adapters ligated together without insert) waste sequencing capacity

**Low Complexity Regions**: Sometimes library preparation produces fragments with unusual sequence composition:
- Poly-A or poly-T stretches from RNA contamination
- Repetitive sequences that may align to multiple genomic locations
- PCR artifacts that create false duplicates

### The Biology of ChIP-seq Fragment Sizes

Understanding the biological basis of ChIP-seq fragment sizes helps guide processing decisions:

**Nucleosome Protection**: In eukaryotic cells, DNA is wrapped around histone proteins forming nucleosomes. Each nucleosome protects approximately 147 base pairs of DNA from nuclease digestion. Even when studying transcription factors like ARF27, the chromatin context means that ChIP-seq fragments often reflect nucleosome-sized protection.

**Cross-linking Range**: Formaldehyde cross-linking can capture protein-DNA interactions within a range of ~200-500 base pairs. This means that even direct transcription factor binding sites may yield fragments of varying sizes depending on:
- Local chromatin structure
- Presence of other proteins in the complex
- Efficiency of sonication during chromatin fragmentation

**Optimal Fragment Size for Analysis**: For transcription factor ChIP-seq like ARF27, fragments between 100-300 base pairs provide the best resolution for identifying precise binding sites while maintaining sufficient signal for detection.

## Trimmomatic: The Workhorse of Read Processing

Trimmomatic is specifically designed for Illumina sequencing data and provides sophisticated algorithms for removing technical artifacts while preserving biological signal.

### Understanding Trimmomatic's Approach

Trimmomatic uses a multi-step approach to clean reads:

1. **Adapter Detection and Removal**: Uses seed-based alignment to find adapter sequences
2. **Quality-based Trimming**: Removes low-quality regions using sliding window approaches
3. **Length Filtering**: Discards reads that become too short after trimming

This order is crucial because:
- Adapters must be removed before quality assessment (they may have different quality patterns)
- Quality trimming after adapter removal is more accurate
- Length filtering at the end ensures you don't waste time processing reads that will ultimately be discarded

### Setting Up Read Processing

```bash
# Ensure we're in the correct environment and directory
cd chipseq_analysis
conda activate chipseq

# Create directory for processed reads
mkdir -p processed_data
```

### Step 1: Check Your Data Type

Before trimming, let's see what type of data you have:

```bash
# Look at your raw data files to understand the naming pattern
ls raw_data/

# Check if you have paired-end data (files ending with _R1 and _R2)
ls raw_data/*_R1.fastq.gz 2>/dev/null || echo "No paired-end R1 files found"
ls raw_data/*_R2.fastq.gz 2>/dev/null || echo "No paired-end R2 files found"
```

### Step 2: Single-End Read Trimming

If your data is single-end (most common for ChIP-seq), process each file individually:

```bash
# Example: Process the first sample (replace with your actual filename)
trimmomatic SE -threads 4 \
    raw_data/sample1.fastq.gz \
    processed_data/sample1_trimmed.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33
```

Repeat this command for each of your samples, changing the filenames:

```bash
# Process sample 2
trimmomatic SE -threads 4 \
    raw_data/sample2.fastq.gz \
    processed_data/sample2_trimmed.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33

# Process sample 3
trimmomatic SE -threads 4 \
    raw_data/sample3.fastq.gz \
    processed_data/sample3_trimmed.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33

# Continue for all your samples...
```

### Step 3: Paired-End Read Trimming (if applicable)

If your data is paired-end, you need to process both R1 and R2 files together:

```bash
# Example: Process paired-end sample 1
trimmomatic PE -threads 4 \
    raw_data/sample1_R1.fastq.gz \
    raw_data/sample1_R2.fastq.gz \
    processed_data/sample1_R1_paired.fastq.gz \
    processed_data/sample1_R1_unpaired.fastq.gz \
    processed_data/sample1_R2_paired.fastq.gz \
    processed_data/sample1_R2_unpaired.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33
```

Repeat for each paired-end sample:

```bash
# Process paired-end sample 2
trimmomatic PE -threads 4 \
    raw_data/sample2_R1.fastq.gz \
    raw_data/sample2_R2.fastq.gz \
    processed_data/sample2_R1_paired.fastq.gz \
    processed_data/sample2_R1_unpaired.fastq.gz \
    processed_data/sample2_R2_paired.fastq.gz \
    processed_data/sample2_R2_unpaired.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33
```

### Step 4: Check What Files Were Created

After trimming, verify your output:

```bash
# See what trimmed files were created
ls -lh processed_data/

# Count how many files you should have
echo "Raw files:"
ls raw_data/*.fastq.gz | wc -l

echo "Processed files:"
ls processed_data/*.fastq.gz | wc -l
```

## Detailed Parameter Explanation

Let's break down each Trimmomatic parameter and understand why it's chosen:

### ILLUMINACLIP Parameters
```
ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10
```

**Adapter File**: `TruSeq3-SE.fa` contains sequences of Illumina TruSeq adapters
- This file includes the most common adapter sequences used in Illumina sequencing
- Using the wrong adapter file will fail to remove contamination
- The conda installation provides standard adapter files for most protocols

**Seed Mismatches (2)**: Allows up to 2 mismatches in the initial 16-base seed match
- **Why 2?**: Balances sensitivity (finding real adapters with sequencing errors) vs. specificity (avoiding false matches)
- Too low (0-1): May miss adapters with sequencing errors
- Too high (3+): May incorrectly identify genomic sequences as adapters

**Palindrome Clip Threshold (30)**: For paired-end data, how well forward and reverse reads must match
- **Why 30?**: Requires strong evidence of adapter read-through
- Higher values: More stringent, fewer false positives
- Lower values: More sensitive but may clip valid sequences

**Simple Clip Threshold (10)**: Minimum score for adapter match in single-end mode
- **Why 10?**: Removes clear adapter contamination while preserving ambiguous cases
- This threshold prevents over-aggressive trimming of sequences that might be genuinely genomic

### Quality-Based Trimming Parameters

**LEADING:3** - Remove low-quality bases from the beginning of reads
- **Why 3?**: Phred score 3 = ~50% accuracy, clearly unreliable
- Removes bases with quality scores below 3 from the 5' end
- These often result from sequencing startup artifacts

**TRAILING:3** - Remove low-quality bases from the end of reads
- **Why 3?**: Same reasoning as LEADING
- More important than LEADING because quality typically degrades toward 3' end
- Prevents poor-quality tails from causing alignment problems

**SLIDINGWINDOW:4:20** - Advanced quality trimming using a sliding window approach
- **Window Size (4)**: Examines 4 consecutive bases at a time
- **Quality Threshold (20)**: Average quality in the window must be ≥20 (99% accuracy)

**Why This Approach Works**:
- More sophisticated than simple "cut at first bad base"
- Allows temporary quality dips that recover
- Removes sustained regions of poor quality
- Window size of 4 provides good balance between sensitivity and stability

**MINLEN:36** - Discard reads shorter than 36 bases after trimming
- **Why 36?**: Minimum length for reliable genome alignment
- Shorter reads become ambiguous (may align to multiple locations)
- 36 bases provides sufficient uniqueness for most genomic regions
- Balances read retention vs. alignment accuracy

### Threading and Performance

**-threads 4**: Uses 4 CPU cores for parallel processing
- **Why 4?**: Good balance for most systems without overwhelming other processes
- Can be increased to 8 or more on high-end systems
- Monitor CPU usage during processing to optimize

**-phred33**: Specifies quality score encoding
- Modern Illumina data uses Phred+33 encoding
- Older data might use Phred+64 (now rare)
- Using wrong encoding leads to incorrect quality interpretation

### Running the Trimming Process

The trimming commands above will process your files one by one. Here's what you'll see:

```bash
# Example output from trimmomatic
TrimmomaticSE: Started with arguments:
 -threads 4 raw_data/sample1.fastq.gz processed_data/sample1_trimmed.fastq.gz 
 ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36 -phred33

Input Reads: 25000000 Surviving: 23456789 (93.82%) Dropped: 1543211 (6.18%)
TrimmomaticSE: Completed successfully
```

**What to Expect During Processing**:
- Processing time: 10-30 minutes per million reads
- Memory usage: Typically <2GB RAM
- Disk space: Processed files are usually 70-90% the size of originals

### Monitoring Progress

While trimming is running, you can check progress:

```bash
# Check how many files have been processed
ls processed_data/*.fastq.gz | wc -l

# Check file sizes (processed files should be smaller than originals)
ls -lh processed_data/

# Monitor system resources
top
```

## Understanding Trimming Statistics

Trimmomatic provides detailed statistics about what it removed:

```
Input Reads: 25000000 Surviving: 23456789 (93.82%) Dropped: 1543211 (6.18%)
TrimmomaticSE: Started with arguments:...
Input Reads: 25000000
Surviving Reads: 23456789
Dropped Reads: 1543211
Surviving Read Length: 65.4 (average)
```

**Key Metrics to Evaluate**:

**Survival Rate**: Percentage of reads retained after processing
- **>95%**: Excellent, high-quality data with minimal contamination
- **90-95%**: Good, some quality issues but well-handled
- **80-90%**: Acceptable, moderate quality problems corrected
- **<80%**: Concerning, may indicate systematic quality problems

**Average Length Reduction**: How much shorter reads became
- **0-5 bases**: Minimal trimming needed, high-quality data
- **5-15 bases**: Moderate trimming, typical for most datasets
- **15-25 bases**: Significant trimming, quality issues present
- **>25 bases**: Extensive trimming, consider data quality carefully

## Post-Trimming Quality Assessment

After trimming, verify that the processing improved data quality:

### Step 1: Run FastQC on Trimmed Data

```bash
# Create directory for trimmed data quality reports
mkdir -p quality_control/trimmed_fastqc

# Run FastQC on all trimmed files
fastqc processed_data/*_trimmed.fastq.gz -o quality_control/trimmed_fastqc/ -t 4
```

For paired-end data, also check the paired files:

```bash
# For paired-end data, check the paired files
fastqc processed_data/*_paired.fastq.gz -o quality_control/trimmed_fastqc/ -t 4
```

### Step 2: Generate Quality Reports

```bash
# Generate report for trimmed data only
multiqc quality_control/trimmed_fastqc/ -o quality_control/ --filename trimmed_data_report

# Create before/after comparison report
multiqc quality_control/raw_fastqc/ quality_control/trimmed_fastqc/ \
    -o quality_control/ --filename before_after_trimming_comparison
```

### Step 3: Compare Before and After

Open these HTML reports in your web browser:
- `quality_control/raw_data_report.html` (from Document 2)
- `quality_control/trimmed_data_report.html` (new)
- `quality_control/before_after_trimming_comparison.html` (comparison)

### Evaluating Trimming Success

Compare the before and after FastQC reports:

**Quality Score Improvements**:
- Per base quality should show improvement, especially at 3' end
- More reads should have high average quality scores
- Quality distribution should shift toward higher values

**Adapter Content Reduction**:
- Adapter contamination should be reduced to <1%
- Any remaining adapter content should be at very low levels
- Complete elimination isn't always possible but dramatic reduction is expected

**Sequence Length Distribution**:
- Should show the expected reduction in average length
- Very short reads should be eliminated
- Distribution should be tighter around the mean

## When Standard Trimming Isn't Enough

### More Aggressive Trimming (if needed)

If your post-trimming quality is still not good enough, try more stringent parameters:

```bash
# Example with more aggressive quality trimming
trimmomatic SE -threads 4 \
    raw_data/sample1.fastq.gz \
    processed_data/sample1_high_quality.fastq.gz \
    ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 \
    LEADING:10 \
    TRAILING:10 \
    SLIDINGWINDOW:4:25 \
    MINLEN:50 \
    -phred33
```

**Changes made**:
- `LEADING:10` and `TRAILING:10`: Remove more bases from ends
- `SLIDINGWINDOW:4:25`: Higher quality threshold (25 instead of 20)
- `MINLEN:50`: Longer minimum length (50 instead of 36)

### Custom Adapter Sequences (if needed)

If you know your specific adapter sequences, create a custom file:

```bash
# Create custom adapter file
cat > custom_adapters.fa << 'EOF'
>TruSeq_Adapter_1
AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
>TruSeq_Adapter_2
AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
EOF

# Use custom adapters in trimming
trimmomatic SE -threads 4 \
    raw_data/sample1.fastq.gz \
    processed_data/sample1_custom_trimmed.fastq.gz \
    ILLUMINACLIP:custom_adapters.fa:2:30:10 \
    LEADING:3 \
    TRAILING:3 \
    SLIDINGWINDOW:4:20 \
    MINLEN:36 \
    -phred33
```

## Optimizing Parameters for ChIP-seq

### Fragment Size Considerations

For transcription factor ChIP-seq like ARF27:

**Minimum Length**: Should preserve the core binding region
- 36 bp minimum is standard but 50 bp may be better for large genomes
- Consider that alignment algorithms need sufficient unique sequence
- Shorter reads increase mapping ambiguity in repetitive regions

**Quality Threshold**: Balance between quality and read retention
- SLIDINGWINDOW:4:20 is standard but consider 4:25 for critical experiments
- Higher thresholds reduce false positive alignments
- Lower thresholds preserve more data but may increase noise

### Maize-Specific Considerations

**Repetitive Content**: Maize's high repetitive content means:
- Higher quality requirements for confident alignment
- Consider slightly more stringent trimming
- Minimum length of 50 bp may be beneficial

**Genome Size**: Large genome means:
- More potential alignment locations for short reads
- Higher quality requirements for unique mapping
- Processing time is longer but trimming remains fast

## Documentation and Quality Tracking

### Create a Simple Processing Summary

```bash
# Create a summary of what you did
cat > quality_control/trimming_summary.txt << 'EOF'
TRIMMING PROCESS SUMMARY
========================

Date: 
Trimmomatic parameters used:
- ILLUMINACLIP: TruSeq3-SE.fa:2:30:10
- LEADING: 3
- TRAILING: 3  
- SLIDINGWINDOW: 4:20
- MINLEN: 36

Sample Processing Results:
EOF
```

### Check How Many Reads Survived

```bash
# Count reads in original files
echo "Original read counts:"
for file in raw_data/*.fastq.gz; do
    sample=$(basename $file .fastq.gz)
    count=$(zcat $file | wc -l)
    reads=$((count / 4))
    echo "$sample: $reads reads"
done

# Count reads in trimmed files  
echo "Trimmed read counts:"
for file in processed_data/*_trimmed.fastq.gz; do
    sample=$(basename $file _trimmed.fastq.gz)
    count=$(zcat $file | wc -l)
    reads=$((count / 4))
    echo "$sample: $reads reads"
done
```

### Calculate Survival Rates

For each sample, you can calculate what percentage of reads survived:

```bash
# Example for one sample (replace sample1 with your actual sample name)
original=$(zcat raw_data/sample1.fastq.gz | wc -l | awk '{print $1/4}')
trimmed=$(zcat processed_data/sample1_trimmed.fastq.gz | wc -l | awk '{print $1/4}')
survival=$(echo "scale=2; $trimmed * 100 / $original" | bc)
echo "Sample1 survival rate: $survival%"
```

### Quality Metrics to Track

Monitor these metrics throughout processing:

1. **Read Survival Rate**: Percentage of reads retained
2. **Average Quality Improvement**: Change in mean quality scores
3. **Length Distribution**: How read lengths changed
4. **Adapter Removal Efficiency**: Reduction in adapter contamination
5. **Processing Time**: For optimization of future runs

## Preparing for Alignment

After successful read processing, you should have:

**Clean Read Files**: High-quality, adapter-free sequences ready for alignment
**Quality Documentation**: Records of what was removed and why
**Processing Statistics**: Data on read survival and quality improvement
**Validated Improvement**: Confirmation that trimming achieved its goals

The processed reads are now ready for genome alignment, where we'll map each read to its location in the maize genome to begin identifying where ARF27 binds.

### Pre-Alignment Checklist

Before proceeding to alignment:

✓ **Quality Verification**: Post-trimming FastQC shows improved quality
✓ **Read Counts**: Sufficient reads remain for robust analysis (>20M for ChIP samples)
✓ **File Organization**: Processed files are properly named and organized
✓ **Documentation**: Processing parameters and results are recorded
✓ **Adapter Removal**: Adapter contamination reduced to <1%

This careful read processing establishes the foundation for accurate genome alignment and ultimately reliable identification of ARF27 binding sites.
