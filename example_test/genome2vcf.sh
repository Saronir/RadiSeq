#!/bin/bash
# Usage: ./genome2vcf.sh reference.fa mutated.fa
# Output: variants.vcf.gz

set -euo pipefail

# Input arguments
REFERENCE=$1
MUTATED=$2

# Output files
ALN_BAM="aln.sorted.bam"
VCF="variants.vcf.gz"
CHUNK_DIR="chunks_tmp"
CHUNK_SIZE=10000000   # 10 Mbp per chunk

# Tools required: minimap2, samtools, bcftools, seqkit
for tool in minimap2 samtools bcftools seqkit; do
    if ! command -v $tool &> /dev/null; then
        echo "❌ Error: $tool not found in PATH. Please install it."
        exit 1
    fi
done

# Function to get genome size
genome_size() {
    awk '/^>/ {if (seqlen){print seqlen}; seqlen=0; next} {seqlen += length($0)} END {print seqlen}' "$1" \
    | awk '{sum+=$1} END {print sum}'
}

# Estimate genome size of mutated genome
GENOME_SIZE=$(genome_size "$MUTATED")
echo ">>> Mutated genome size ≈ $GENOME_SIZE bp"

# If genome is huge (>100 Mbp), use chunking
if [ "$GENOME_SIZE" -gt 100000000 ]; then
    echo ">>> Genome is large, using chunking (chunks of $CHUNK_SIZE bp)..."
    mkdir -p $CHUNK_DIR
    seqkit split -s $CHUNK_SIZE "$MUTATED" -O $CHUNK_DIR

    echo ">>> Aligning chunks..."
    for CHUNK in $CHUNK_DIR/*.fa; do
        minimap2 -a "$REFERENCE" "$CHUNK" | samtools view -bS -
    done | samtools sort -o $ALN_BAM

    rm -rf $CHUNK_DIR
else
    echo ">>> Genome is small, aligning directly..."
    minimap2 -a "$REFERENCE" "$MUTATED" \
        | samtools view -bS - \
        | samtools sort -o $ALN_BAM
fi

echo ">>> Indexing BAM"
samtools index $ALN_BAM

echo ">>> Calling variants"
bcftools mpileup -f "$REFERENCE" $ALN_BAM \
  | bcftools call -mv -Oz -o $VCF

echo ">>> Indexing VCF"
bcftools index $VCF

echo "✅ Done! Variants written to: $VCF"
