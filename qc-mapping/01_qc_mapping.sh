#!/bin/bash

# Project: Different sex determination systems in two Phoxinus species
# Pipeline: FastQC → fastp → BWA Mapping → SAMtools → Picard → Qualimap
# Cluster: Leibniz Institute for the Analysis of Biodiversity Change (LIB)
# Scheduler: Sun Grid Engine
# Author: Temitope Oriowo, 2025
# Requirements: fastqc, fastp, bwa-mem2, samtools, picard, qualimap, java
# Notes: Assumes input FASTQ files are in subdirectories under ../Fastq/

#==================== SETUP ====================
REFERENCE="fPhoPho.hap1.fa"
THREADS=16

FASTQ_DIR="../Fastq"
TRIM_DIR="../Trimmed_Quality"
QC_DIR="../Quality_Check"
MAP_DIR="../Mapping"
TMP_DIR="../tmp"
PICARD_JAR="../Tools/picard.jar"
QUALIMAP_OUT="${MAP_DIR}/Qualimap_results"

# Trimming/filter thresholds
MIN_LENGTH=25
TRIM_FRONT=5      # trim from 5' end
TRIM_TAIL=5       # trim from 3' end
MIN_QUAL=20

mkdir -p "$QC_DIR" "$TRIM_DIR" "$MAP_DIR" "$TMP_DIR" "$QUALIMAP_OUT"

#==================== LOAD MODULES ====================
module load fastqc/0.12.1
module load fastp/0.23.2
module load bwa2/2.2.1
module load samtools/1.19.2
module load picard/3.2.0
module load java/jre1.8.0_231
module load qualimap/2.3

#==================== STEP 1: FASTQC ====================
echo "Starting FastQC..."
for subdir in "$FASTQ_DIR"/*/; do
    sample_dir=$(basename "$subdir")
    mkdir -p "$QC_DIR/$sample_dir"
    cd "$subdir" || continue

    for file in *_2.fq.gz; do
        sample="${file/_2.fq.gz/}"
        if [[ -f "${sample}_1.fq.gz" ]]; then
            fastqc "${sample}_1.fq.gz" "${sample}_2.fq.gz" --outdir="$QC_DIR/$sample_dir"
            echo "FastQC complete for ${sample}"
        fi
    done
done

#==================== STEP 2: TRIMMING WITH fastp ====================
echo "Starting fastp trimming..."
for subdir in "$FASTQ_DIR"/*/; do
    sample_dir=$(basename "$subdir")
    mkdir -p "$TRIM_DIR/$sample_dir"
    cd "$subdir" || continue

    for file in *_1.fq.gz; do
        sample="${file/_1.fq.gz/}"
        paired_file="${sample}_2.fq.gz"

        if [[ -f "$paired_file" ]]; then
            fastp \
                --in1 "$file" \
                --in2 "$paired_file" \
                --out1 "${TRIM_DIR}/${sample_dir}/${sample}_1_trimmed.fq" \
                --out2 "${TRIM_DIR}/${sample_dir}/${sample}_2_trimmed.fq" \
                --length_required "$MIN_LENGTH" \
                --cut_front --cut_front_window_size 1 --cut_front_mean_quality "$MIN_QUAL" \
                --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality "$MIN_QUAL" \
                --trim_front1 "$TRIM_FRONT" --trim_tail1 "$TRIM_TAIL" \
                --trim_front2 "$TRIM_FRONT" --trim_tail2 "$TRIM_TAIL" \
                --qualified_quality_phred "$MIN_QUAL" \
                --thread "$THREADS" \
                --html "${TRIM_DIR}/${sample_dir}/${sample}_fastp.html" \
                --json "${TRIM_DIR}/${sample_dir}/${sample}_fastp.json"

            echo "fastp trimming complete for ${sample}"
        fi
    done
done

#==================== STEP 3: MAPPING ====================
echo "Indexing reference and starting BWA mapping..."
bwa-mem2 index "$REFERENCE" || exit 1

cd "$TRIM_DIR" || exit 1

for sample_dir in */; do
    cd "$sample_dir" || continue

    for file in *_1_trimmed.fq; do
        sample="${file/_1_trimmed.fq/}"
        RG="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA"

        fq1="${sample}_1_trimmed.fq"
        fq2="${sample}_2_trimmed.fq"
        sam="${MAP_DIR}/${sample}.sam"
        bam="${MAP_DIR}/${sample}.bam"
        sorted_bam="${MAP_DIR}/${sample}.sorted.bam"
        dedup_bam="${MAP_DIR}/${sample}.sorted.rm.bam"
        dup_metrics="${MAP_DIR}/${sample}_duplication_metrics.txt"

        bwa-mem2 mem -M -R "$RG" -t "$THREADS" "$REFERENCE" "$fq1" "$fq2" > "$sam" || continue

        samtools view -b -q 30 "$sam" > "$bam" && rm "$sam"

        ulimit -n 2048
        samtools sort -T "$TMP_DIR" -o "$sorted_bam" "$bam" && rm "$bam"

        java -jar "$PICARD_JAR" MarkDuplicates \
            I="$sorted_bam" \
            O="$dedup_bam" \
            M="$dup_metrics" \
            REMOVE_DUPLICATES=true \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 || continue

        unset DISPLAY
        qualimap bamqc -bam "$dedup_bam" -outdir "$QUALIMAP_OUT/${sample}" \
            --java-mem-size=12G -nt 1 -outfile "${sample}_qualimap_stats.pdf"

        echo "$(date '+%Y-%m-%d %H:%M:%S') Finished ${sample}" >> "${MAP_DIR}/mapping.log"
    done

    cd ..
done

echo "Pipeline completed successfully!"
