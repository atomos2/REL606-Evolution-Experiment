set -e

echo "Running FastQC..."
cd ..

output_dir=results/FastQC
mkdir -p "$output_dir"

fastqc data/*.fastq.gz -o "$output_dir"

echo "Unzipping FastQC .zip files..."
cd "$output_dir"

for file in *.zip; do
    unzip -o "$file"
done

summary_dir=../../docs
mkdir -p "$summary_dir"

echo "Creating summary for all files..."
cat */summary.txt > "$summary_dir/summaries.txt"

echo "FastQC processing complete."
  