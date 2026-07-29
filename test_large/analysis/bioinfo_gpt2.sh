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

# Structural variant options.
# Manta is the main recommended SV caller.
RUN_MANTA="${RUN_MANTA:-1}"

# Delly is optional. Set RUN_DELLY=1 if installed.
RUN_DELLY="${RUN_DELLY:-1}"

# Use the Bowtie2 BAM for SV calling by default.
# Set USE_BWA_FOR_SV=1 to create a separate BWA-MEM BAM for Manta/Delly.
USE_BWA_FOR_SV="${USE_BWA_FOR_SV:-0}"

mkdir -p "${OUTDIR}"
mkdir -p "$(dirname "${INDEX}")"

# -------------------------------------------------------------------------
# Output files
# -------------------------------------------------------------------------

BAM="${OUTDIR}/${SAMPLE}.bowtie2.sorted.bam"

RAW_BCF="${OUTDIR}/${SAMPLE}.raw.bcf"
NORM_BCF="${OUTDIR}/${SAMPLE}.normalized.bcf"
FILTERED_VCF="${OUTDIR}/${SAMPLE}.filtered.vcf.gz"
PASS_VCF="${OUTDIR}/${SAMPLE}.pass.vcf.gz"

BWA_BAM="${OUTDIR}/${SAMPLE}.bwa.sorted.bam"

MANTA_DIR="${OUTDIR}/manta"
MANTA_DIPLOID_VCF="${MANTA_DIR}/results/variants/diploidSV.vcf.gz"
MANTA_CANDIDATE_VCF="${MANTA_DIR}/results/variants/candidateSV.vcf.gz"
MANTA_PASS_VCF="${OUTDIR}/${SAMPLE}.manta.pass.vcf.gz"
MANTA_DIPLOID_SUMMARY="${OUTDIR}/${SAMPLE}.manta_diploid_sv_summary.tsv"
MANTA_CANDIDATE_SUMMARY="${OUTDIR}/${SAMPLE}.manta_candidate_sv_summary.tsv"

DELLY_BCF="${OUTDIR}/${SAMPLE}.delly.bcf"
DELLY_VCF="${OUTDIR}/${SAMPLE}.delly.vcf.gz"
DELLY_SUMMARY="${OUTDIR}/${SAMPLE}.delly_sv_summary.tsv"

# -------------------------------------------------------------------------
# Check required programs
# -------------------------------------------------------------------------

for PROGRAM in bowtie2 bowtie2-build samtools bcftools; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} was not found in PATH." >&2
        exit 1
    fi
done

if [[ "${RUN_MANTA}" == "1" ]]; then
    if ! command -v configManta.py >/dev/null 2>&1; then
        echo "ERROR: configManta.py was not found in PATH." >&2
        echo "Install Manta or set RUN_MANTA=0 to skip Manta SV calling." >&2
        exit 1
    fi
fi

if [[ "${USE_BWA_FOR_SV}" == "1" ]]; then
    if ! command -v bwa >/dev/null 2>&1; then
        echo "ERROR: bwa was not found in PATH." >&2
        echo "Install bwa or set USE_BWA_FOR_SV=0." >&2
        exit 1
    fi
fi

if [[ "${RUN_DELLY}" == "1" ]]; then
    if ! command -v delly >/dev/null 2>&1; then
        echo "WARNING: delly was not found in PATH. Delly SV calling will be skipped." >&2
        RUN_DELLY=0
    fi
fi

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
    echo

    bowtie2 --version | head -n 1 || true
    samtools --version | head -n 1 || true
    bcftools --version | head -n 1 || true

    if command -v bwa >/dev/null 2>&1; then
        echo
        { bwa 2>&1 | head -n 3; } || true
    fi

    if command -v configManta.py >/dev/null 2>&1; then
        echo
        echo "Manta configManta.py found at: $(command -v configManta.py)"
    fi

    if command -v delly >/dev/null 2>&1; then
        echo
        { delly --version 2>&1; } || true
    fi
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
# Bowtie2 alignment and coordinate sorting
# -------------------------------------------------------------------------

echo "Aligning ${SAMPLE} with Bowtie2..."

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
samtools quickcheck -v "${BAM}"

# -------------------------------------------------------------------------
# Optional BWA-MEM alignment branch for SV calling
# -------------------------------------------------------------------------

SV_BAM="${BAM}"

if [[ "${USE_BWA_FOR_SV}" == "1" ]]; then
    echo "Building BWA-MEM alignment branch for SV calling..."

    if [[ ! -s "${REF}.bwt" ]]; then
        echo "Building BWA index..."
        bwa index "${REF}"
    fi

    bwa mem \
        -t "${THREADS}" \
        -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:RadiSeq\tPL:ILLUMINA" \
        "${REF}" \
        "${R1}" \
        "${R2}" \
    | samtools sort \
        -@ "${THREADS}" \
        -m 1G \
        -O BAM \
        -o "${BWA_BAM}" \
        -

    samtools index -@ "${THREADS}" "${BWA_BAM}"
    samtools quickcheck -v "${BWA_BAM}"

    SV_BAM="${BWA_BAM}"
fi

echo "BAM used for SNP/indel calling: ${BAM}"
echo "BAM used for SV calling:        ${SV_BAM}"

# -------------------------------------------------------------------------
# Alignment quality control
# -------------------------------------------------------------------------

echo "Generating alignment QC..."

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
# SNP and short-indel calling with bcftools
# -------------------------------------------------------------------------

echo "Calling SNPs and short indels with bcftools..."

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
# Soft filtering for SNPs and short indels
# -------------------------------------------------------------------------

bcftools filter \
    -s LowQual \
    -e 'QUAL<20 || FORMAT/DP<5 || FORMAT/AD[0:1]<2' \
    -Oz \
    -o "${FILTERED_VCF}" \
    "${NORM_BCF}"

bcftools index -f "${FILTERED_VCF}"

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

# -------------------------------------------------------------------------
# Structural variant calling with Manta
# -------------------------------------------------------------------------

if [[ "${RUN_MANTA}" == "1" ]]; then
    echo "Calling structural variants with Manta..."

    rm -rf "${MANTA_DIR}"

    configManta.py \
        --bam "${SV_BAM}" \
        --referenceFasta "${REF}" \
        --runDir "${MANTA_DIR}"

    "${MANTA_DIR}/runWorkflow.py" \
        -m local \
        -j "${THREADS}"

    if [[ ! -s "${MANTA_DIPLOID_VCF}" ]]; then
        echo "ERROR: Manta did not produce expected file: ${MANTA_DIPLOID_VCF}" >&2
        exit 1
    fi

    if [[ ! -s "${MANTA_CANDIDATE_VCF}" ]]; then
        echo "ERROR: Manta did not produce expected file: ${MANTA_CANDIDATE_VCF}" >&2
        exit 1
    fi

    bcftools index -f "${MANTA_DIPLOID_VCF}" || true
    bcftools index -f "${MANTA_CANDIDATE_VCF}" || true

    # PASS-only diploid Manta calls.
    bcftools view \
        -f PASS \
        -Oz \
        -o "${MANTA_PASS_VCF}" \
        "${MANTA_DIPLOID_VCF}"

    bcftools index -f "${MANTA_PASS_VCF}"

    # Diploid SV summary.
    bcftools query \
        -u \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/MATEID\t%INFO/EVENT[\t%GT\t%PR\t%SR]\n' \
        "${MANTA_DIPLOID_VCF}" \
        > "${MANTA_DIPLOID_SUMMARY}"

    # Candidate SV summary. Useful for simulator benchmarking because it is more sensitive.
    bcftools query \
        -u \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/MATEID\t%INFO/EVENT[\t%GT\t%PR\t%SR]\n' \
        "${MANTA_CANDIDATE_VCF}" \
        > "${MANTA_CANDIDATE_SUMMARY}"

    echo "Manta diploid SV VCF:      ${MANTA_DIPLOID_VCF}"
    echo "Manta candidate SV VCF:    ${MANTA_CANDIDATE_VCF}"
    echo "Manta PASS SV VCF:         ${MANTA_PASS_VCF}"
    echo "Manta diploid summary:     ${MANTA_DIPLOID_SUMMARY}"
    echo "Manta candidate summary:   ${MANTA_CANDIDATE_SUMMARY}"
fi

# -------------------------------------------------------------------------
# Optional structural variant calling with Delly
# -------------------------------------------------------------------------

if [[ "${RUN_DELLY}" == "1" ]]; then
    echo "Calling structural variants with Delly..."

    delly sr \
        -g "${REF}" \
        -o "${DELLY_BCF}" \
        "${SV_BAM}"

    bcftools view \
        -Oz \
        -o "${DELLY_VCF}" \
        "${DELLY_BCF}"

    bcftools index -f "${DELLY_VCF}"

    bcftools query \
        -u \
        -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/END\t%INFO/SVLEN\t%INFO/CHR2\t%INFO/POS2[\t%GT]\n' \
        "${DELLY_VCF}" \
        > "${DELLY_SUMMARY}"

    echo "Delly SV VCF:      ${DELLY_VCF}"
    echo "Delly SV summary:  ${DELLY_SUMMARY}"
fi

# -------------------------------------------------------------------------
# Convenience SV type summaries
# -------------------------------------------------------------------------

if [[ "${RUN_MANTA}" == "1" ]]; then
    echo "Generating Manta SV type count summary..."

    {
        echo "SVTYPE count from Manta candidateSV.vcf.gz"
        bcftools query \
            -u \
            -f '%INFO/SVTYPE\n' \
            "${MANTA_CANDIDATE_VCF}" \
        | sort \
        | uniq -c \
        | sort -nr
    } > "${OUTDIR}/${SAMPLE}.manta_candidate_svtype_counts.txt"

    {
        echo "SVTYPE count from Manta diploidSV.vcf.gz"
        bcftools query \
            -u \
            -f '%INFO/SVTYPE\n' \
            "${MANTA_DIPLOID_VCF}" \
        | sort \
        | uniq -c \
        | sort -nr
    } > "${OUTDIR}/${SAMPLE}.manta_diploid_svtype_counts.txt"
fi

if [[ "${RUN_DELLY}" == "1" ]]; then
    echo "Generating Delly SV type count summary..."

    {
        echo "SVTYPE count from Delly VCF"
        bcftools query \
            -u \
            -f '%INFO/SVTYPE\n' \
            "${DELLY_VCF}" \
        | sort \
        | uniq -c \
        | sort -nr
    } > "${OUTDIR}/${SAMPLE}.delly_svtype_counts.txt"
fi

# -------------------------------------------------------------------------
# Final output summary
# -------------------------------------------------------------------------

echo
echo "Pipeline completed successfully."
echo
echo "Main outputs:"
echo "  Bowtie2 BAM:              ${BAM}"
echo "  SNP/short-indel VCF:      ${FILTERED_VCF}"
echo "  SNP/short-indel PASS VCF: ${PASS_VCF}"
echo

if [[ "${USE_BWA_FOR_SV}" == "1" ]]; then
    echo "  BWA SV BAM:               ${BWA_BAM}"
fi

if [[ "${RUN_MANTA}" == "1" ]]; then
    echo "  Manta diploid SV VCF:     ${MANTA_DIPLOID_VCF}"
    echo "  Manta candidate SV VCF:   ${MANTA_CANDIDATE_VCF}"
    echo "  Manta PASS SV VCF:        ${MANTA_PASS_VCF}"
    echo "  Manta candidate summary:  ${MANTA_CANDIDATE_SUMMARY}"
fi

if [[ "${RUN_DELLY}" == "1" ]]; then
    echo "  Delly SV VCF:             ${DELLY_VCF}"
    echo "  Delly SV summary:         ${DELLY_SUMMARY}"
fi

echo
echo "Useful notes:"
echo "  DEL = deletion"
echo "  INV = inversion"
echo "  BND or TRA = translocation/breakend-style event"
echo "  candidateSV.vcf.gz is often better for simulator truth comparison than PASS-only calls."
