# Long-Term *E. coli* Evolution Experiment – Variant Calling and Hypermutation Analysis

## Overview

This project investigates **hypermutability** in *Escherichia coli* populations from the **Lenski Long-Term Evolution Experiment (LTEE)**. The main objective is to detect genomic variants across three evolutionary timepoints (5,000, 15,000, and 50,000 generations) and identify mutations in **DNA mismatch repair genes**—particularly *mutL*—that may underlie the emergence of hypermutability.

## Background

Richard Lenski's LTEE is a landmark experiment studying **evolutionary adaptation** by propagating 12 replicate populations of *E. coli* strain REL606 in a glucose-limited minimal medium for over 40,000 generations. One of these populations (Ara-3) evolved the ability to metabolize citrate aerobically (Cit+), and some lineages acquired **hypermutator phenotypes**, often linked to mutations in **mismatch repair (MMR)** genes (*mutL*, *mutS*, *mutH*).

In this analysis we used sequencing data from three timepoints in the Ara-3 population:

* **SRR2589044** — 5,000 generations
* **SRR2584863** — 15,000 generations
* **SRR2584866** — 50,000 generations

## Tools Required

* `FastQC` – for quality control
* `Trimmomatic` – for adapter trimming
* `Bowtie2` – for read alignment
* `Samtools` – for format conversion and indexing
* `Bcftools` – for variant calling

## Step 1: Downloading Data

```bash
cd ../data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/004/SRR2589044/SRR2589044_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/003/SRR2584863/SRR2584863_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR258/006/SRR2584866/SRR2584866_2.fastq.gz
```

## Step 2: Quality Check (FastQC)

```bash
mkdir -p results/FastQC
fastqc data/*.fastq.gz -o results/FastQC
cd results/FastQC
for file in *.zip; do unzip -o "$file"; done
cat */summary.txt > ../../docs/summaries.txt
```

To view the summaries:

```bash
less ../../docs/summaries.txt
firefox *.html
```

## Step 3: Trimming Adapters (Trimmomatic)

```bash
mkdir -p results/trimmed results/orphaned
for fastq in data/*_1.fastq.gz; do
  filename=$(basename "$fastq" _1.fastq.gz)
  TrimmomaticPE \
    data/${filename}_1.fastq.gz data/${filename}_2.fastq.gz \
    results/trimmed/${filename}_1.trimmed.fastq.gz results/orphaned/${filename}_1.untrimmed.fastq.gz \
    results/trimmed/${filename}_2.trimmed.fastq.gz results/orphaned/${filename}_2.untrimmed.fastq.gz \
    SLIDINGWINDOW:4:20 MINLEN:25 \
    ILLUMINACLIP:/usr/share/trimmomatic/NexteraPE-PE.fa:2:40:15
done
```

## Step 4: Variant Calling Workflow

```bash
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
```

## Step 5: Variant Counts Per Generation

```bash
grep -v "#" results/vcf/SRR2589044_final_variants.vcf | wc -l   # 29 variants
grep -v "#" results/vcf/SRR2584863_final_variants.vcf | wc -l   # 32 variants
grep -v "#" results/vcf/SRR2584866_final_variants.vcf | wc -l   # 829 variants
```

## Step 6: Extract *mutL* Variants

(*Gene location: NC\_012967.1:4375567–4377414*)

```bash
VCF_DIR="results/vcf"
CHROM="NC_012967.1"
START=4375567
END=4377414
OUTDIR="results/mutLvar"
mkdir -p "$OUTDIR"

for vcf in ${VCF_DIR}/*_final_variants.vcf; do
  sample=$(basename "$vcf" _final_variants.vcf)
  if [[ ! -f "$vcf.gz" ]]; then bgzip -c "$vcf" > "$vcf.gz"; fi
  if [[ ! -f "$vcf.gz.tbi" ]]; then tabix -p vcf "$vcf.gz"; fi
  bcftools view -r ${CHROM}:${START}-${END} "$vcf.gz" -o "${OUTDIR}/${sample}_mutL.vcf"
done
```

## Final Result & Conclusion

* **SRR2589044 (5,000 gen)** → 29 variants
* **SRR2584863 (15,000 gen)** → 32 variants
* **SRR2584866 (50,000 gen)** → 829 variants

Only **SRR2584866** (50,000 generations) contains a **mutation in the *mutL* gene** at **position 4,377,265** (*A → G*).
This strongly suggests that **hypermutability** emerged due to the disruption of mismatch repair via the *mutL* mutation.

**Conclusion:** Mutation in *mutL* correlates with the rise in mutation rate in the 50,000-generation sample.

## Credits

* Data Source: [ENA Project PRJNA188723](https://www.ebi.ac.uk/ena/browser/view/PRJNA188723)
* Tutorial Reference: [GTK Teaching - NGS Intro](https://gtk-teaching.github.io/NGS-intro/)

## Author

**Iftikhar Alam**
MSc Bioinformatics, Pondicherry University
