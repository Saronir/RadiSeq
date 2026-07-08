#!/usr/bin/env bash

set -Eeuo pipefail

trap 'echo "ERROR: Pipeline failed at line ${LINENO}" >&2' ERR

# -------------------------------------------------------------------------
# Configuration
# -------------------------------------------------------------------------

THREADS="${THREADS:-8}"
SAMPLE="${SAMPLE:-CELL_0}"

# IMPORTANT:
# REF must be the ORIGINAL, UNMUTATED reference genome.
REF="${REF:-../radiseq_humanlike_reference.fa}"

# Bowtie 2 index prefix built from the same REF file.
INDEX="${INDEX:-./index/synthetic_genome}"

R1="${R1:-../output/${SAMPLE}_R1.fastq}"
R2="${R2:-../output/${SAMPLE}_R2.fastq}"

OUTDIR="${OUTDIR:-results/${SAMPLE}}"

# Set according to the simulated fragment-size distribution.
MAX_INSERT="${MAX_INSERT:-1000}"

mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${INDEX}")"

BAM="${OUTDIR}/${SAMPLE}.sorted.bam"
RAW_BCF="${OUTDIR}/${SAMPLE}.raw.bcf"
NORM_BCF="${OUTDIR}/${SAMPLE}.normalized.bcf"
FILTERED_VCF="${OUTDIR}/${SAMPLE}.filtered.vcf.gz"
PASS_VCF="${OUTDIR}/${SAMPLE}.pass.vcf.gz"

# -------------------------------------------------------------------------
# Check required programs
# -------------------------------------------------------------------------

for PROGRAM in bowtie2 bowtie2-build samtools bcftools; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} was not found in PATH." >&2
        exit 1
    fi
done

# -------------------------------------------------------------------------
# Check input files
# -------------------------------------------------------------------------

for FILE in "${REF}" "${R1}" "${R2}"; do
    if [[ ! -s "${FILE}" ]]; then
        echo "ERROR: Input file is missing or empty: ${FILE}" >&2
        exit 1
    fi
done

# -------------------------------------------------------------------------
# Record software versions
# -------------------------------------------------------------------------

{
    echo "Pipeline run: $(date --iso-8601=seconds)"
    bowtie2 --version | head -n 1
    samtools --version | head -n 1
    bcftools --version | head -n 1
} > "${OUTDIR}/software_versions.txt"

# -------------------------------------------------------------------------
# Build reference indexes if required
# -------------------------------------------------------------------------

if [[ ! -s "${REF}.fai" ]]; then
    echo "Building FASTA index..."
    samtools faidx "${REF}"
fi

if [[ ! -e "${INDEX}.1.bt2" && ! -e "${INDEX}.1.bt2l" ]]; then
    echo "Building Bowtie 2 index..."
    bowtie2-build "${REF}" "${INDEX}"
fi

# -------------------------------------------------------------------------
# Alignment and coordinate sorting
# -------------------------------------------------------------------------

echo "Aligning ${SAMPLE}..."

bowtie2 \
    --very-sensitive-local \
    --fr \
    -X "${MAX_INSERT}" \
    -p "${THREADS}" \
    -x "${INDEX}" \
    -1 "${R1}" \
    -2 "${R2}" \
    --rg-id "${SAMPLE}" \
    --rg "SM:${SAMPLE}" \
    --rg "LB:RadiSeq" \
    --rg "PL:ILLUMINA" \
    2> "${OUTDIR}/${SAMPLE}.bowtie2.log" \
| samtools sort \
    -@ "${THREADS}" \
    -m 1G \
    -O BAM \
    -o "${BAM}" \
    -

samtools index -@ "${THREADS}" "${BAM}"

# Verify that the BAM has a valid header and EOF block.
samtools quickcheck -v "${BAM}"

# -------------------------------------------------------------------------
# Alignment quality control
# -------------------------------------------------------------------------

samtools flagstat \
    -@ "${THREADS}" \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.flagstat.txt"

samtools stats \
    -@ "${THREADS}" \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.samtools_stats.txt"

samtools coverage \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.coverage.tsv"

samtools idxstats \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.idxstats.tsv"

# Optional per-base depth file. This can be very large for a full genome.
samtools depth \
    -a \
    "${BAM}" \
    > "${OUTDIR}/${SAMPLE}.depth.tsv"

# -------------------------------------------------------------------------
# SNP and short-indel calling
# -------------------------------------------------------------------------

echo "Calling SNPs and short indels..."

bcftools mpileup \
    -Ou \
    -f "${REF}" \
    -q 20 \
    -Q 20 \
    -d 10000 \
    -a FORMAT/AD,FORMAT/DP \
    "${BAM}" \
| bcftools call \
    -m \
    -v \
    --ploidy 2 \
    -Ob \
    -o "${RAW_BCF}"

bcftools index -f "${RAW_BCF}"

# -------------------------------------------------------------------------
# Variant normalization
# -------------------------------------------------------------------------

bcftools norm \
    -f "${REF}" \
    -m -any \
    -Ob \
    -o "${NORM_BCF}" \
    "${RAW_BCF}"

bcftools index -f "${NORM_BCF}"

# -------------------------------------------------------------------------
# Soft filtering
# -------------------------------------------------------------------------

# Starting thresholds:
#   QUAL >= 20
#   sample depth >= 5
#   alternate allele supported by at least 2 reads
#
# Failed calls are retained and labelled "LowQual".

bcftools filter \
    -s LowQual \
    -e 'QUAL<20 || FORMAT/DP<5 || FORMAT/AD[0:1]<2' \
    -Oz \
    -o "${FILTERED_VCF}" \
    "${NORM_BCF}"

bcftools index -f "${FILTERED_VCF}"

# Generate a second file containing only PASS variants.
bcftools view \
    -f PASS \
    -Oz \
    -o "${PASS_VCF}" \
    "${FILTERED_VCF}"

bcftools index -f "${PASS_VCF}"

# -------------------------------------------------------------------------
# Variant quality control
# -------------------------------------------------------------------------

bcftools stats \
    "${FILTERED_VCF}" \
    > "${OUTDIR}/${SAMPLE}.variant_stats.txt"

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL[\t%GT\t%DP\t%AD]\n' \
    "${FILTERED_VCF}" \
    > "${OUTDIR}/${SAMPLE}.variant_summary.tsv"

echo "Pipeline completed successfully."
echo "Aligned BAM:       ${BAM}"
echo "All filtered calls: ${FILTERED_VCF}"
echo "PASS calls:         ${PASS_VCF}"