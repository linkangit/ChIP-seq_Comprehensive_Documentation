# Step 1: Read Quality Control for ChIP-seq Data

## Objective
Assess the quality of raw ChIP-seq FASTQ files to identify potential issues before mapping and ensure high-quality downstream analysis.

## Why Quality Control is Critical

Quality control is the foundation of any successful ChIP-seq analysis. Poor quality data will lead to:
- Low mapping rates
- False positive peaks
- Reduced sensitivity to detect true binding sites
- Unreliable biological conclusions

ARF27 ChIP-seq data quality is especially important because transcription factors typically have:
- Lower signal-to-noise ratios compared to histone marks
- Fewer binding sites (hundreds to thousands vs. millions for histones)
- More susceptible to technical artifacts

## What We're Looking For

### Good Quality Indicators:
- Per-base quality scores >20 (preferably >30)
- Even GC content distribution
- Low adapter contamination
- Minimal overrepresented sequences
- Appropriate read length distribution

### Warning Signs:
- Quality drops significantly toward 3' end
- High adapter content
- Unusual GC content bias
- High duplication rates (though some is expected in ChIP-seq)
- Overrepresented sequences (could be contamination)

## Step-by-Step Instructions

### 1. Set Up Your Working Directory

```bash
# Create organized directory structure
mkdir -p chipseq_arf27_analysis
cd chipseq_arf27_analysis

# Create subdirectories for each analysis step
mkdir -p 01_raw_data
mkdir -p 02_qc_reports
mkdir -p 03_trimmed_data
mkdir -p 04_mapped_data
mkdir -p 05_processed_bams
mkdir -p 06_peaks
mkdir -p 07_annotations
mkdir -p 08_motifs
mkdir -p 09_visualization
mkdir -p logs

# Move your FASTQ files to the raw_data directory
# Example file names (adjust according to your data):
# ARF27_ChIP_rep1_R1.fastq.gz, ARF27_ChIP_rep1_R2.fastq.gz
# ARF27_ChIP_rep2_R1.fastq.gz, ARF27_ChIP_rep2_R2.fastq.gz
# Input_control_rep1_R1.fastq.gz, Input_control_rep1_R2.fastq.gz
# Input_control_rep2_R1.fastq.gz, Input_control_rep2_R2.fastq.gz
```

### 2. Install FastQC (if not already installed)

```bash
# Using conda (recommended)
conda install -c bioconda fastqc

# Or download from: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
# Extract and add to PATH
```

### 3. Run FastQC on All FASTQ Files

```bash
# Navigate to raw data directory
cd 01_raw_data

# Run FastQC on all FASTQ files
# This will process each file individually
fastqc ARF27_ChIP_rep1_R1.fastq.gz
fastqc ARF27_ChIP_rep1_R2.fastq.gz
fastqc ARF27_ChIP_rep2_R1.fastq.gz
fastqc ARF27_ChIP_rep2_R2.fastq.gz
fastqc Input_control_rep1_R1.fastq.gz
fastqc Input_control_rep1_R2.fastq.gz
fastqc Input_control_rep2_R1.fastq.gz
fastqc Input_control_rep2_R2.fastq.gz

# Alternative: Run on all files at once
# fastqc *.fastq.gz

# Move HTML reports to QC directory
mv *.html ../02_qc_reports/
mv *.zip ../02_qc_reports/
```

**Command Explanation:**
- `fastqc`: The quality control tool
- `*.fastq.gz`: Processes all gzipped FASTQ files
- Output: HTML reports (human-readable) and ZIP files (machine-readable)

### 4. Examine FastQC Reports

Open the HTML reports in a web browser. For each sample, examine these key sections:

#### A. Basic Statistics
Check:
- Total sequences (should be similar across replicates)
- Sequence length (typically 50-150 bp)
- %GC content (should be ~42% for maize, may vary slightly for ChIP samples)

#### B. Per Base Sequence Quality
- **Good**: Quality scores >30 across most of the read
- **Acceptable**: Quality >20, with some decline toward 3' end
- **Problematic**: Quality <20 for significant portions

#### C. Per Sequence Quality Scores
- Most reads should have average quality >25
- Very few reads with quality <15

#### D. Per Base Sequence Content
- Should be relatively even across positions
- Some bias at the beginning is normal due to random priming
- Extreme bias may indicate adapter contamination

#### E. Adapter Content
- **Good**: <5% adapter content throughout the read
- **Needs trimming**: >10% adapter content, especially toward 3' end

#### F. Overrepresented Sequences
- Small amounts are normal in ChIP-seq
- High levels may indicate PCR bias or contamination
- Check if sequences match known adapters or contaminants

## Expected Output

You should see HTML reports for each FASTQ file with graphical summaries. Here's what to expect:

### Typical ChIP-seq Characteristics:
- **Moderate duplication**: ChIP-seq naturally has more duplicates than RNA-seq
- **GC bias**: May differ from genome average due to protein binding preferences
- **Read quality**: Should be high (>Q20) for most of the read

## Quality Assessment Decision Tree

### If Quality is Good (Most metrics pass):
```bash
# Proceed directly to mapping (Step 2)
echo "Quality looks good - proceeding to mapping"
```

### If Adapter Contamination is Present:
```bash
# Plan to trim adapters (we'll cover this in mapping step)
echo "Need to trim adapters before mapping"
```

### If Quality Drops Significantly:
```bash
# Plan to trim low-quality bases
echo "Need quality trimming"
```

## Common Issues and Troubleshooting

### Issue 1: Very High Duplication Rates (>80%)
**Possible Causes:**
- Over-amplification during PCR
- Low complexity library
- Very successful IP (high enrichment)

**Solutions:**
- Check if this is consistent across replicates
- If only one replicate affected, consider excluding it
- Proceed with duplicate removal (Step 4)

**Command to check duplication:**
```bash
# Count total reads
zcat ARF27_ChIP_rep1_R1.fastq.gz | wc -l | awk '{print $1/4}'

# This gives total read count - compare with unique reads after duplicate removal
```

### Issue 2: Adapter Contamination
**Identification:**
- FastQC shows high adapter content
- Overrepresented sequences match adapter sequences

**Solutions:**
- Note which adapters are present
- Plan adapter trimming during mapping

### Issue 3: Poor Quality Toward 3' End
**Normal for:** Older sequencing platforms
**Solutions:**
- Plan quality trimming
- Consider if read length after trimming will be sufficient (>30bp recommended)

### Issue 4: Unusual GC Content
**For ARF27 ChIP-seq:**
- ARF proteins bind AT-rich auxin response elements
- Slight AT bias compared to genome average is expected
- Major deviations may indicate contamination

### Issue 5: Very Low Read Counts
**Minimum recommendations:**
- ChIP samples: 20-50 million reads per replicate
- Input controls: 15-30 million reads per replicate
- For transcription factors like ARF27: aim for higher end of range

## Quality Control Checklist

Before proceeding to Step 2, verify:

- [ ] All FASTQ files processed successfully
- [ ] HTML reports generated and examined
- [ ] Overall quality scores are acceptable (>Q20 average)
- [ ] Noted any adapter contamination for trimming
- [ ] Documented any quality issues in analysis log
- [ ] Read counts are sufficient for analysis

## Analysis Log Template

Keep detailed notes of your findings:

```bash
# Create analysis log
cat > ../logs/step1_qc_log.txt << EOF
Step 1: Quality Control Results
Date: $(date)

Sample Quality Summary:
- ARF27_ChIP_rep1: [Good/Needs trimming/Poor] - Notes: 
- ARF27_ChIP_rep2: [Good/Needs trimming/Poor] - Notes:
- Input_control_rep1: [Good/Needs trimming/Poor] - Notes:
- Input_control_rep2: [Good/Needs trimming/Poor] - Notes:

Issues Identified:
- Adapter contamination: [Yes/No]
- Quality trimming needed: [Yes/No]
- Any samples to exclude: [None/List samples]

Decisions for Next Step:
- Proceed to mapping: [Yes/No]
- Trimming parameters needed: [List if applicable]

EOF
```

## Next Steps

Based on your QC results:
1. **If quality is good**: Proceed to Step 2 (Read Mapping)
2. **If trimming needed**: Note parameters for trimming during mapping
3. **If major issues**: Consider re-sequencing problematic samples

## Key Takeaways

- Quality control is not optional - it guides all downstream decisions
- ChIP-seq has different quality expectations than RNA-seq
- Document everything for reproducible analysis
- When in doubt, err on the side of caution with quality thresholds

Ready to move on to Step 2: Read Mapping? This is where we'll align your quality-controlled reads to the maize B73-v4 reference genome.
