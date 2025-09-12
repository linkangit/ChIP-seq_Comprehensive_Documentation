# Document 2: Quality Control and Data Assessment

## The Critical Importance of Quality Control

Quality control is the foundation of reliable ChIP-seq analysis. Poor quality data will produce unreliable results regardless of how sophisticated your downstream analysis becomes. Think of QC as the diagnostic phase of your analysis - you're examining your data to understand its strengths, weaknesses, and potential problems before investing time in complex analyses.

In ChIP-seq, quality issues can arise at multiple stages:
- **Library preparation**: Adapter contamination, fragment size issues, low complexity
- **Sequencing**: Base calling errors, uneven coverage, failed chemistry
- **Sample preparation**: Poor cross-linking, inefficient immunoprecipitation, degraded DNA

Identifying these issues early allows you to either correct them computationally or recognize when experimental problems require new samples.

## Understanding Your Raw Data

### What FASTQ Files Contain

Your raw sequencing data comes in FASTQ format, which contains four lines per sequence read:

```
@HWI-D00119:50:H7AP8ADXX:1:1101:1234:2000 1:N:0:ATCACG
GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
+
!''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
```

**Line 1**: Read identifier containing sequencer information, run details, and barcode
**Line 2**: The actual DNA sequence (A, T, G, C, N)
**Line 3**: Separator (just a "+")
**Line 4**: Quality scores for each base (encoded as ASCII characters)

### Quality Score Encoding

Quality scores represent the confidence that each base was called correctly. They're encoded using ASCII characters where:
- Each character represents a quality score (Phred score)
- Phred score 20 = 99% accuracy (1 in 100 chance of error)
- Phred score 30 = 99.9% accuracy (1 in 1000 chance of error)
- Phred score 40 = 99.99% accuracy (1 in 10,000 chance of error)

**Why This Matters**: Low quality bases contribute noise to your analysis. Bases with quality scores below 20 are essentially unreliable and should be trimmed or filtered out.

## Running Initial Quality Assessment

### Setting Up Quality Control

```bash
# Go to your analysis folder
cd chipseq_analysis

# Create a folder for quality results
mkdir quality_control
```

### FastQC Analysis - Your First Look at Data Quality

```bash
# Run quality check on your files
fastqc *.fastq.gz -o quality_control/
```

**What this does**: FastQC checks your sequencing files and creates quality reports for each one.

### Understanding FastQC Output

FastQC generates both HTML reports (human-readable) and text files (machine-readable). Let's examine the key analyses:

#### 1. Basic Statistics
```
Filename: sample1.fastq.gz
File type: Conventional base calls
Encoding: Sanger / Illumina 1.9
Total Sequences: 25,000,000
Sequences flagged as poor quality: 0
Sequence length: 75
%GC: 42
```

**What to Look For**:
- **Total Sequences**: Should match expected sequencing depth (20-60 million for ChIP-seq)
- **Sequence length**: Should be consistent with your sequencing protocol
- **%GC**: Should be reasonable for your organism (maize genome is ~46% GC)

#### 2. Per Base Sequence Quality
This is arguably the most important plot. It shows quality scores across all positions in your reads.

**Green Zone (Quality 28+)**: Excellent quality, very reliable bases
**Orange Zone (Quality 20-28)**: Reasonable quality, acceptable for most analyses  
**Red Zone (Quality <20)**: Poor quality, should be trimmed or filtered

**Normal Patterns**:
- Quality typically starts high and declines toward the 3' end
- Single-end reads often show quality drop after position 50-60
- Paired-end reads may show quality dips in the middle (due to sequencing chemistry)

**Warning Signs**:
- Quality dropping below 20 for large portions of reads
- Unusual spikes or dips in quality
- Very poor quality from the beginning (suggests sequencing problems)

#### 3. Per Sequence Quality Scores
Shows the distribution of average quality scores across all reads.

**Good Pattern**: Most reads should have average quality scores above 25
**Warning Signs**: 
- Large numbers of reads with low average quality
- Bimodal distributions (suggesting mixed quality populations)

#### 4. Per Base Sequence Content
Shows the proportion of each nucleotide (A, T, G, C) at each position.

**Expected Pattern**: 
- Lines should be roughly horizontal and parallel
- %A should approximately equal %T
- %G should approximately equal %C
- Some variation at the beginning is normal (due to random priming)

**Warning Signs**:
- Strong bias toward specific nucleotides
- Dramatic changes in composition along read length
- Extreme GC bias (could indicate contamination)

#### 5. Sequence Duplication Levels
Shows what percentage of your reads are duplicates.

**ChIP-seq Specifics**: Unlike RNA-seq, ChIP-seq naturally has some duplication due to:
- PCR amplification during library preparation
- True biological enrichment (same DNA fragments pulled down multiple times)
- Limited complexity of enriched regions

**Acceptable Levels**: 
- <30% duplication: Excellent
- 30-50% duplication: Good for ChIP-seq
- 50-70% duplication: Acceptable but may reduce peak resolution
- Over 70% duplication: Concerning, may indicate over-amplification

#### 6. Adapter Content
Shows contamination with sequencing adapters.

**Why This Matters**: Adapters are artificial sequences added during library preparation. If not properly removed, they:
- Reduce the amount of useful genomic sequence
- Can cause alignment problems
- Create false signals in downstream analysis

**Normal Pattern**: Should be 0% or very low across all positions
**Warning Signs**: Adapter content increasing toward 3' end of reads

### Creating Comprehensive Quality Reports with MultiQC

Individual FastQC reports are detailed but can be overwhelming when analyzing multiple samples. MultiQC aggregates results into a single, comparative report.

```bash
# Create a combined report for all samples
multiqc quality_control/ -o quality_control/
```

**What MultiQC Provides**:
- Side-by-side comparison of all samples
- Summary statistics across the entire dataset
- Identification of outlier samples
- Interactive plots for detailed exploration

**Key Sections to Review**:

1. **General Statistics Table**: Overview of all samples showing read counts, duplication levels, GC content
2. **Sequence Quality Histograms**: Compare quality distributions across samples
3. **Per Sequence GC Content**: Identify samples with unusual GC bias
4. **Sequence Duplication Levels**: Spot samples with excessive duplication

## Interpreting Quality Results

### Sample-Level Assessment

For each sample, ask these questions:

**Is the sequencing depth adequate?**
- ChIP samples: 20-60 million reads
- Input controls: 10-40 million reads (can be lower than ChIP)

**Is the sequence quality sufficient?**
- Most bases should have quality scores >20
- If quality drops dramatically, you'll need aggressive trimming

**Are there technical artifacts?**
- High adapter content requires trimming
- Unusual sequence composition might indicate contamination
- Extreme duplication levels suggest over-amplification

### Dataset-Level Assessment

**Are samples comparable?**
- Similar read counts across biological replicates
- Consistent quality metrics between related samples
- No obvious outliers that might skew analysis

**Are ChIP and input samples properly paired?**
- Input samples should have lower duplication than ChIP (less enrichment)
- Similar read depths help with normalization
- Quality metrics should be comparable

### Red Flags That Require Action

**Immediate Concerns** (may require new sequencing):
- Very low read counts (<10 million)
- Extremely poor quality (most bases <20)
- High contamination levels
- Complete absence of expected sequences

**Correctable Issues** (can fix with trimming/filtering):
- Moderate quality decline toward 3' end
- Adapter contamination
- Some low-quality reads mixed with good ones

**ChIP-seq Specific Issues**:
- Input samples with higher duplication than ChIP samples (suggests swapped labels)
- Dramatically different read counts between replicates
- Unusual GC content patterns that don't match expected organism

## Quality Control Metrics for ChIP-seq Success

### Expected Quality Patterns in ChIP-seq

**ChIP Samples**:
- Moderate duplication levels (30-60%) due to enrichment
- Relatively normal GC content (close to genome average)
- Good overall sequence quality
- Consistent metrics between biological replicates

**Input/Control Samples**:
- Lower duplication levels than ChIP samples
- GC content closer to genome background
- Similar quality metrics to ChIP samples
- May have slightly fewer reads (input is less "interesting" to sequence deeply)

### Cross-Sample Comparisons

Use MultiQC to identify:

**Outlier Samples**: Those with dramatically different metrics from others
**Batch Effects**: Systematic differences between sequencing runs
**Label Swaps**: Input samples that look like ChIP samples (or vice versa)

## Documentation and Decision Making

### Creating a Quality Assessment Report

Document your findings systematically:

```bash
# Create a simple summary file
echo "Quality Assessment Summary" > quality_summary.txt
echo "Date: $(date)" >> quality_summary.txt
echo "Total samples: [fill in number]" >> quality_summary.txt
echo "Samples needing trimming: [list here]" >> quality_summary.txt
echo "Samples with issues: [list here]" >> quality_summary.txt
```

### Decision Points Based on Quality Assessment

**Proceed with standard analysis** if:
- Most samples have >20 million reads
- Quality scores generally >20
- Duplication levels reasonable for ChIP-seq
- No major technical artifacts

**Proceed with modified analysis** if:
- Some quality issues that can be corrected with trimming
- Moderate technical artifacts
- Some samples with borderline quality

**Consider re-sequencing** if:
- Very low read counts across samples
- Severe quality problems that can't be corrected
- Major contamination issues
- Critical samples completely failed

### Quality Control as an Iterative Process

Quality control doesn't end here. You'll continue monitoring data quality throughout the analysis:

**After trimming**: Verify that quality improvement was achieved
**After alignment**: Check mapping rates and alignment quality
**After peak calling**: Assess peak quality and reproducibility between replicates

Each step provides new quality metrics that inform your interpretation of results.

## Setting Quality Standards

### Minimum Acceptable Standards

**Read Quality**: 
- At least 80% of bases with quality ≥20
- Mean read quality ≥25

**Read Count**:
- ChIP samples: ≥20 million reads
- Input samples: ≥10 million reads

**Technical Quality**:
- Adapter contamination <5%
- Duplication levels <80%

### Optimal Standards

**Read Quality**:
- At least 95% of bases with quality ≥20
- Mean read quality ≥30

**Read Count**:
- ChIP samples: 40-60 million reads
- Input samples: 20-40 million reads

**Technical Quality**:
- Adapter contamination <1%
- Duplication levels 30-60%

## Common Quality Issues and Solutions

### Low Quality Scores
**Cause**: Sequencing chemistry problems, old reagents, overloaded flowcells
**Solution**: Aggressive quality trimming, potentially re-sequencing if severe

### High Adapter Content
**Cause**: Insert sizes shorter than read length, incomplete adapter removal
**Solution**: Adapter trimming with appropriate parameters

### Unusual GC Content
**Cause**: Contamination, PCR bias, species mismatch
**Solution**: Investigate source, potentially exclude contaminated samples

### Extreme Duplication
**Cause**: Over-amplification, low library complexity, very strong enrichment
**Solution**: Evaluate during peak calling, may need to adjust duplicate handling

### Inconsistent Quality Between samples
**Cause**: Different sequencing runs, batch effects, sample degradation
**Solution**: Note for downstream analysis, may need batch correction

## Preparing for the Next Step

Based on your quality assessment, you now know:

1. **Which samples are suitable for analysis**
2. **What trimming parameters to use**
3. **Whether any samples need special handling**
4. **What quality issues to monitor in downstream steps**

This information directly informs the read processing and trimming step, where you'll clean up the data to maximize the quality of your downstream analysis.

The time invested in thorough quality control pays dividends throughout the rest of your analysis by ensuring you're working with the best possible data and understanding its limitations.
