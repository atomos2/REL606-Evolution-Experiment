
set -e
cd ..
mkdir -p results/trimmed
mkdir -p results/orphaned

echo "Trimmomatic getting started. . ."

for fastq in data/*_1.fastq.gz; do

                filename=$(basename "$fastq" _1.fastq.gz)
                echo "Running trimmomatic for ${filename}..."

                TrimmomaticPE \
                        data/${filename}_1.fastq.gz data/${filename}_2.fastq.gz \
                        results/trimmed/${filename}_1.trimmed.fastq.gz results/orphaned/${filename}_1.untrimmed.fastq.gz \
                        results/trimmed/${filename}_2.trimmed.fastq.gz results/orphaned/${filename}_2.untrimmed.fastq.gz \
                        SLIDINGWINDOW:4:20 MINLEN:25 \
                        ILLUMINACLIP:/usr/share/trimmomatic/NexteraPE-PE.fa:2:40:15
done


echo "Trimming done for all samples"