#!/bin/bash
#standard bioinformatics pipeline from Samy's work

bowtie2 -x ./index -1 ../output/CELL_2_R1.fastq -2 ../output/CELL_2_R2.fastq -S outsam.sam

samtools view -b -S ./outsam.sam -o outbam.bam

samtools sort -o outsorted.bam ./outbam.bam

samtools index ./outsorted.bam

bcftools mpileup -Ou -f ../contig2.fa ./outsorted.bam --threads 1 | bcftools call -m -Ov | bcftools norm -m -any | bcftools filter -e 'QUAL<10 || DP<10' > outfiltered.vcf
