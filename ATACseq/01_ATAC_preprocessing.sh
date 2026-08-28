#!/bin/bash
# ATAC-seq Data Preprocessing Pipeline
# =============================================================================
# This script performs a complete ATAC-seq analysis workflow including:
#   1. Quality control (FastQC, MultiQC)
#   2. Adapter trimming (fastp)
#   3. Alignment to reference genome (Bowtie2)
#   4. BAM file processing and filtering (samtools, Picard)
#   5. Coverage file generation (deepTools)
#   6. Peak calling (MACS2)
#
# USAGE:
#   1. Prepare your raw FASTQ files in the ./rawdata/ directory
#   2. Create a file named 'name.txt' containing sample names (one per line)
#      Sample naming format: sample_R1_*.fastq.gz and sample_R2_*.fastq.gz
#   3. Modify the configuration variables below (paths, genome version, etc.)
#   4. Run: bash atac_seq_pipeline.sh
#
# REQUIREMENTS:
#   - Conda environment with: multiqc, macs2, deepTools, subread, FastQC, fastp, Bowtie2, samtools, Picard, bedtools (paths configured below)
#   - Reference genome Bowtie2 index
#   - Blacklist regions BED file for your genome
#   - Gene annotation BED file
#
# OUTPUT STRUCTURE:
#   ./analysis_results/
#   ├── 01_pre_qc/          - Pre-trimming FastQC reports
#   ├── 02_trimmed/         - Trimmed FASTQ files
#   ├── 03_post_qc/         - Post-trimming FastQC reports
#   ├── 04_alignment/       - BAM files and alignment statistics
#   ├── 05_coverage/        - BigWig coverage files and heatmaps
#   └── 06_peaks/           - MACS2 peak calling results
#
# NOTE: Read shifting (Tn5 offset correction) is NOT performed in this pipeline.
# =============================================================================

set -euo pipefail

# =============================================================================
# CONFIGURATION SECTION - Modify these variables for your analysis
# =============================================================================

# Computing Resources
THREADS=32                          # Number of CPU threads to use for parallel processing

# Input/Output Directories
RAW_DATA_DIR="./rawdata"          # Directory containing input FASTQ files
OUTPUT_DIR="./analysis_results"   # Main output directory for all results

# Analysis Parameters
MAX_frag_len=1000                   # Maximum fragment length for Bowtie2 alignment
                                    # Fragments longer than this will be discarded

# Genome Size for MACS2 peak calling
# Effective genome sizes for uniquely mappable regions (150bp reads):
#   - GRCh38 (human): 2862010428
#   - mm10 (mouse):   2494787038
GENOME_SIZE="2494787038"          # Change to 2862010428 for human GRCh38

# =============================================================================
# REFERENCE FILE PATHS - Update these to match your system
# =============================================================================
declare -A REF_PATHS=(
    # Bowtie2 index prefix (without .1.bt2, .2.bt2, etc. extensions)
    ["bowtie2_index"]="${HOME}/datapool/rawdata/ref/Bowtie2/mm10_forATAC/mm10"

    # ENCODE blacklist regions BED file (regions to exclude from analysis)
    ["blacklist"]="${HOME}/datapool/rawdata/ref/blacklist/mm10-blacklist.v2.bed"

    # Gene annotation BED file for TSS heatmap generation
    # Format: chr start end gene_name strand
    ["genes_bed"]="${HOME}/datapool/rawdata/ref/bed/mm10_Refseq.bed"
)

# =============================================================================
# SOFTWARE PATHS - Update these to match your installation locations
# =============================================================================
declare -A SOFTWARE_PATHS=(
    ["fastqc"]="${HOME}/datapool/soft/FastQC/fastqc"
    ["fastp"]="${HOME}/datapool/soft/fastp"
    ["bowtie2"]="${HOME}/datapool/soft/bowtie2-2.5.4/bowtie2"
    ["samtools"]="${HOME}/datapool/soft/samtools-1.20/samtools"
    ["picard"]="${HOME}/datapool/soft/picard.jar"
    ["bedtools"]="${HOME}/datapool/soft/bedtools2/bin/bedtools"
)

# =============================================================================
# OUTPUT DIRECTORY STRUCTURE
# =============================================================================
declare -A DIRS=(
    ["pre_qc"]="${OUTPUT_DIR}/01_pre_qc"         # Pre-trimming quality control
    ["trimmed"]="${OUTPUT_DIR}/02_trimmed"       # Trimmed FASTQ files
    ["post_qc"]="${OUTPUT_DIR}/03_post_qc"       # Post-trimming quality control
    ["alignment"]="${OUTPUT_DIR}/04_alignment"   # Alignment results
    ["coverage"]="${OUTPUT_DIR}/05_coverage"     # Coverage tracks and heatmaps
    ["peaks"]="${OUTPUT_DIR}/06_peaks"           # Peak calling results
)

# =============================================================================
# MAIN PIPELINE START
# =============================================================================

echo -e "\n\e[34m====== ATAC-seq Analysis Pipeline Started ======\e[0m"

# Check for required sample list file
if [ ! -f name.txt ]; then
  echo "ERROR: name.txt not found!"
  echo "Please create a file named 'name.txt' containing all sample names (without suffix), one per line."
  echo "Example content of name.txt:"
  echo "  sample1"
  echo "  sample2"
  echo "  sample3"
  exit 1
fi

# Create output directories
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${DIRS[@]}"

# =============================================================================
# STEP 1: Quality Control and Adapter Trimming
# =============================================================================
echo -e "\e[32mStep 1: Quality Control and Trimming...\e[0m"

# Pre-trimming QC: Assess raw data quality with FastQC
# Output: HTML reports and ZIP archives for each FASTQ file
"${SOFTWARE_PATHS[fastqc]}" -t "${THREADS}" -o "${DIRS[pre_qc]}" "${RAW_DATA_DIR}"/*.fastq.gz

# Aggregate FastQC reports with MultiQC
multiqc "${DIRS[pre_qc]}"/*.zip -o "${DIRS[pre_qc]}"/

# Adapter trimming with fastp
# This loop processes each sample listed in name.txt
while read sample; do

    # Find all R1 and R2 files for this sample using pattern matching
    r1_pattern="${RAW_DATA_DIR}/${sample}_R1_*.fastq.gz"
    r1_files=( $r1_pattern )

    r2_pattern="${RAW_DATA_DIR}/${sample}_R2_*.fastq.gz"
    r2_files=( $r2_pattern )

    # Validate that R1 and R2 file counts match
    if [ ${#r1_files[@]} -ne ${#r2_files[@]} ]; then
        echo "ERROR: Mismatched number of R1 and R2 files for sample $sample"
        echo "  R1 files found: ${#r1_files[@]}"
        echo "  R2 files found: ${#r2_files[@]}"
        exit 1
    fi

    # Process each pair of R1/R2 files
    for idx in "${!r1_files[@]}"; do

      r1_file="${r1_files[idx]}"
      r2_file="${r2_files[idx]}"

      if [[ "$r1_file" =~ _R1_([0-9]+)\.fastq\.gz$ ]]; then
          suffix="_${BASH_REMATCH[1]}"  # Extract suffix like '_001'
      else
          echo "ERROR: Unable to parse file name for $sample"
          echo "Expected format: ${sample}_R1_001.fastq.gz"
          exit 1
      fi

      "${SOFTWARE_PATHS[fastp]}" \
        --in1 "$r1_file" \
        --in2 "$r2_file" \
        --out1 "${DIRS[trimmed]}/${sample}_R1_${suffix}.clean.fastq.gz" \
        --out2 "${DIRS[trimmed]}/${sample}_R2_${suffix}.clean.fastq.gz" \
        --json "${DIRS[trimmed]}/${sample}_${suffix}_fastp.json" \
        --html "${DIRS[trimmed]}/${sample}_${suffix}_fastp.html" \
        --detect_adapter_for_pe \
        --thread "${THREADS}"

    done

done < name.txt

# Post-trimming QC: Verify trimming effectiveness
echo -e "\e[32mRunning post-trimming quality control...\e[0m"
"${SOFTWARE_PATHS[fastqc]}" -t "${THREADS}" -o "${DIRS[post_qc]}" "${DIRS[trimmed]}"/*.clean.fastq.gz
multiqc "${DIRS[post_qc]}"/*.zip -o "${DIRS[post_qc]}"/

# =============================================================================
# STEP 2: Sequence Alignment and BAM Processing
# =============================================================================
echo -e "\e[32mStep 2: Alignment and Filtering...\e[0m"

while read sample; do

  echo "Processing sample: $sample"

  r1_pattern="${DIRS[trimmed]}/${sample}_R1_*.clean.fastq.gz"
  r1_files=( $r1_pattern )

  r2_pattern="${DIRS[trimmed]}/${sample}_R2_*.clean.fastq.gz"
  r2_files=( $r2_pattern )

  IFS=, r1_list="${r1_files[*]}" r2_list="${r2_files[*]}"

  # Align reads with Bowtie2
  "${SOFTWARE_PATHS[bowtie2]}" --local --very-sensitive --no-mixed --no-discordant \
    -X "${MAX_frag_len}" -p "${THREADS}" \
    -x "${REF_PATHS[bowtie2_index]}" \
    -1 "${r1_list}" \
    -2 "${r2_list}" \
    2> "${DIRS[alignment]}/${sample}.align.summary" \
    | "${SOFTWARE_PATHS[samtools]}" sort -@ "${THREADS}" -O BAM \
    -o "${DIRS[alignment]}/${sample}.sorted.bam"

  "${SOFTWARE_PATHS[samtools]}" index -@ "${THREADS}" "${DIRS[alignment]}/${sample}.sorted.bam" \
    -o "${DIRS[alignment]}/${sample}.sorted.bam.bai"
  "${SOFTWARE_PATHS[samtools]}" flagstat -@ "${THREADS}" "${DIRS[alignment]}/${sample}.sorted.bam" \
    > "${DIRS[alignment]}/${sample}.sorted.flagstat"
  "${SOFTWARE_PATHS[samtools]}" idxstats -@ "${THREADS}" "${DIRS[alignment]}/${sample}.sorted.bam" \
    > "${DIRS[alignment]}/${sample}.sorted.idxstats"

  # Remove mitochondrial DNA reads (chrM)
  "${SOFTWARE_PATHS[samtools]}" view -h -@ "${THREADS}" "${DIRS[alignment]}/${sample}.sorted.bam" \
    | grep -v "chrM" \
    | "${SOFTWARE_PATHS[samtools]}" sort -@ "${THREADS}" -O bam \
    -o "${DIRS[alignment]}/${sample}.rmChrM.bam"

  # Mark PCR duplicates with Picard
  java -jar "${SOFTWARE_PATHS[picard]}" MarkDuplicates \
    QUIET=true CREATE_INDEX=true REMOVE_DUPLICATES=false \
    VALIDATION_STRINGENCY=LENIENT TMP_DIR=. \
    INPUT="${DIRS[alignment]}/${sample}.rmChrM.bam" \
    OUTPUT="${DIRS[alignment]}/${sample}.rmChrM.markDup.bam" \
    METRICS_FILE="${DIRS[alignment]}/${sample}.dup.metrics"

  # Filter BAM file
  "${SOFTWARE_PATHS[samtools]}" view -h -b -@ "${THREADS}" -f 2 -F 1548 -q 30 \
    "${DIRS[alignment]}/${sample}.rmChrM.markDup.bam" \
    | "${SOFTWARE_PATHS[bedtools]}" intersect -nonamecheck -v \
        -abam - \
        -b "${REF_PATHS[blacklist]}" \
    | "${SOFTWARE_PATHS[samtools]}" sort -@ "${THREADS}" \
    -o "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam"

  "${SOFTWARE_PATHS[samtools]}" index -@ "${THREADS}" \
    "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
    -o "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam.bai"
  "${SOFTWARE_PATHS[samtools]}" flagstat -@ "${THREADS}" \
    "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
    > "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.flagstat"
  "${SOFTWARE_PATHS[samtools]}" idxstats -@ "${THREADS}" \
    "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
    > "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.idxstats"

  # Generate fragment length distribution
  "${SOFTWARE_PATHS[samtools]}" view "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
   | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n \
   | sed -e 's/^[ \t]*//' > "${DIRS[alignment]}/${sample}.fragment_length_count.txt"

done < name.txt

# =============================================================================
# STEP 3: Generate Coverage Files
# =============================================================================
echo -e "\e[32mStep 3: Generating Coverage Files...\e[0m"

while read sample; do

  # Create normalized BigWig coverage track
  bamCoverage --numberOfProcessors "${THREADS}" --binSize 10 \
    --normalizeUsing RPKM --effectiveGenomeSize ${GENOME_SIZE} \
    -b "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
    -o "${DIRS[coverage]}/${sample}_RPKM.bw"

done < name.txt

computeMatrix reference-point --referencePoint TSS --skipZeros \
  -p "${THREADS}" -a 2000 -b 2000 \
  -R "${REF_PATHS[genes_bed]}" \
  -S "${DIRS[coverage]}/"*_RPKM.bw \
  -o "${DIRS[coverage]}/matrix_TSS_RPKM.gz" \
  --outFileSortedRegions "${DIRS[coverage]}/regions_TSS_RPKM.bed"

plotHeatmap -m "${DIRS[coverage]}/matrix_TSS_RPKM.gz" \
  -out "${DIRS[coverage]}/Heatmap_TSS_RPKM.png"

# =============================================================================
# STEP 4: Peak Calling with MACS2
# =============================================================================
echo -e "\e[32mStep 4: Peak Calling...\e[0m"

while read sample; do

  echo "Calling peaks for sample: $sample"

  # Call peaks with MACS2
  macs2 callpeak \
    -f BAMPE --keep-dup all --cutoff-analysis \
    -g "${GENOME_SIZE}" \
    -t "${DIRS[alignment]}/${sample}.rmChrM.markDup.filtered.bam" \
    -n "${sample}" --outdir "${DIRS[peaks]}/${sample}" \
    2> "${DIRS[peaks]}/${sample}.macs2.log"

done < name.txt

# =============================================================================
# PIPELINE COMPLETION
# =============================================================================

echo -e "\n\e[34m====== Analysis Complete! Results saved in ${OUTPUT_DIR} ======\e[0m"
echo ""
echo "Next steps:"
echo "  1. Review MultiQC reports for quality metrics"
echo "  2. Check fragment length distributions for nucleosome patterns"
echo "  3. Visualize peaks in genome browser (load .bw and .narrowPeak files)"
echo "  4. Perform differential accessibility analysis if you have replicates"