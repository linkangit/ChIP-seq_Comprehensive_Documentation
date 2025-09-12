# Beginner-Friendly ChIP-seq Quality Control Guide

## Why Quality Control Matters

Think of quality control like checking your ingredients before cooking. You want to make sure your sequencing data is good quality before spending time analyzing it. Poor quality data will give you unreliable results, no matter how fancy your analysis tools are.

## What You're Looking For

In ChIP-seq data, problems can happen at different stages:
- **During sample prep**: DNA might be degraded or contaminated
- **During sequencing**: The machine might make errors reading the DNA
- **During library prep**: Extra sequences (adapters) might get stuck to your DNA

Finding these problems early helps you decide if you can fix them or if you need new samples.

## Understanding Your Data Files

Your sequencing data comes in FASTQ files. Each DNA sequence has 4 lines:

```
@ReadID_12345
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

- **Line 1**: Name/ID of this DNA piece
- **Line 2**: The actual DNA sequence (A, T, G, C)
- **Line 3**: Just a "+" separator
- **Line 4**: Quality scores (how confident the machine was reading each letter)

## Step 1: Set Up Your Workspace

```bash
# Go to your project folder
cd chipseq_analysis

# Activate your software environment
conda activate chipseq

# Create a folder for quality control results
mkdir quality_control
```

**What this does**: Sets up a clean workspace for your analysis.

## Step 2: Run FastQC (Quality Checker)

FastQC is like a health checkup for your data. It looks at many different aspects of data quality.

```bash
# Check quality of all your data files
fastqc raw_data/*.fastq.gz -o quality_control/
```

**What this command means**:
- `fastqc`: The quality checking program
- `raw_data/*.fastq.gz`: Check all FASTQ files in the raw_data folder
- `-o quality_control/`: Put the results in the quality_control folder

**What happens**: FastQC reads through every DNA sequence in your files and creates detailed reports about quality.

## Step 3: Look at Individual Sample Reports

FastQC creates HTML files you can open in your web browser. Each report has several sections:

### Basic Information
```
Total Sequences: 25,000,000
Sequence length: 75
%GC: 42
```

**What to check**:
- **Total Sequences**: Should be 20-60 million for good ChIP-seq
- **%GC**: Should be close to your organism's normal GC content (about 46% for corn)

### Quality Scores (Most Important!)

Look for the "Per Base Sequence Quality" graph:

- **Green area (score 28+)**: Excellent quality ✅
- **Yellow area (score 20-28)**: OK quality ⚠️
- **Red area (score <20)**: Poor quality, needs fixing ❌

**Normal pattern**: Quality usually starts high and drops toward the end of reads.

**Problems to watch for**:
- Quality dropping below 20 early in the reads
- Weird spikes or dips in quality

### Duplicates

ChIP-seq naturally has some duplicate reads (same DNA piece read multiple times). This is normal because:
- You're enriching for specific DNA regions
- PCR creates copies during library prep

**Acceptable levels**:
- Less than 30%: Excellent
- 30-50%: Good 
- 50-70%: OK but not ideal
- Over 70%: Might be a problem

### Adapters

Adapters are artificial DNA pieces added during library prep. You don't want them in your final data.

**Good**: 0% adapter content
**Bad**: Increasing adapter content, especially at the end of reads

## Step 4: Create a Summary Report

Instead of looking at each sample separately, create one combined report:

```bash
# Create a summary report for all samples
multiqc quality_control/ -o quality_control/ --filename summary_report
```

**What this does**: Combines all individual FastQC reports into one easy-to-read summary.

## Step 5: Interpret Your Results

Open the MultiQC HTML report in your browser. Look for:

### Red Flags (Need to Fix or Re-sequence)
- Very few reads (less than 10 million)
- Most quality scores below 20
- High contamination
- Samples that look very different from others

### Yellow Flags (Can Probably Fix)
- Quality dropping at the end of reads (can trim)
- Some adapter contamination (can remove)
- Moderate quality issues

### Green Flags (Good to Go)
- 20+ million reads
- Most quality scores above 20
- Low adapter contamination
- Similar quality between related samples

## Step 6: Make Decisions

Based on your quality check, decide what to do next:

### If Quality is Good
```bash
# Document your findings
echo "Quality assessment complete. All samples passed quality control." > quality_control/decision.txt
echo "Proceeding with standard analysis pipeline." >> quality_control/decision.txt
```

### If Quality Needs Improvement
```bash
# Note what needs to be fixed
echo "Quality assessment complete. Issues found:" > quality_control/decision.txt
echo "- Need to trim low quality bases from ends" >> quality_control/decision.txt
echo "- Need to remove adapter contamination" >> quality_control/decision.txt
echo "Proceeding with trimming step." >> quality_control/decision.txt
```

### If Quality is Too Poor
```bash
# Document serious problems
echo "Quality assessment complete. Major issues found:" > quality_control/decision.txt
echo "- Several samples have very low read counts" >> quality_control/decision.txt
echo "- Poor quality throughout reads" >> quality_control/decision.txt
echo "Recommend re-sequencing before proceeding." >> quality_control/decision.txt
```

## Quick Quality Checklist

For each sample, check:
- ✅ At least 20 million reads?
- ✅ Most bases have quality score above 20?
- ✅ Reasonable amount of duplicates (less than 70%)?
- ✅ Low adapter contamination (less than 5%)?
- ✅ Similar to other samples in the group?

If you can check most of these boxes, your data is probably good enough to continue.

## Common Problems and Simple Fixes

### Problem: Quality drops at end of reads
**Fix**: Trim the low-quality ends (next step in analysis)

### Problem: Adapter contamination
**Fix**: Remove adapters during trimming step

### Problem: Too many duplicates
**Fix**: Remove duplicates during later processing

### Problem: One sample looks very different
**Check**: Make sure sample labels are correct, might be contaminated

## What's Next?

After quality control, you'll know:
1. Which samples are good enough to analyze
2. What problems need to be fixed in the trimming step
3. Whether any samples should be excluded

Remember: Good quality control at the beginning saves time and frustration later!
