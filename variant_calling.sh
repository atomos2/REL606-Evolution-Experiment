mkdir -p results/sam results/bam results/bcf results/vcf
ref_dir=data/reference/rel606.fna
export BOWTIE2_INDEXES=$(pwd)/data/reference
bowtie2-build "$ref_dir" INrel606

for fastq1 in results/trimmed/*_1.trimmed.fastq.gz; do
  filename=$(basename "$fastq1" _1.trimmed.fastq.gz)
  fastq2=results/trimmed/${filename}_2.trimmed.fastq.gz
  sam=results/sam/${filename}.sam
  bam=results/bam/${filename}.bam
  bam_sorted=results/bam/${filename}_sorted.bam
  bcf=results/bcf/${filename}_raw.bcf
  vcf=results/vcf/${filename}_variants.vcf
  final_vcf=results/vcf/${filename}_final_variants.vcf

  bowtie2 -x INrel606 -1 "$fastq1" -2 "$fastq2" -S "$sam" --very-fast -p 4
  samtools view -S -b "$sam" > "$bam"
  samtools sort -o "$bam_sorted" "$bam"
  samtools index "$bam_sorted"
  bcftools mpileup -O b -o "$bcf" -f "$ref_dir" "$bam_sorted"
  bcftools call --ploidy 1 -m -v -o "$vcf" "$bcf"
  vcfutils.pl varFilter "$vcf" > "$final_vcf"
done
