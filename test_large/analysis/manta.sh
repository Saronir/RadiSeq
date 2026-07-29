SAMPLE=CELL_0
REF=../radiseq_humanlike_reference.fa
BAM=results/${SAMPLE}/${SAMPLE}.bwa.sorted.bam
MANTA_DIR=results/${SAMPLE}/manta

rm -rf "${MANTA_DIR}"

configManta.py \
  --bam "${BAM}" \
  --referenceFasta "${REF}" \
  --runDir "${MANTA_DIR}"

"${MANTA_DIR}/runWorkflow.py" \
  -m local \
  -j 8
