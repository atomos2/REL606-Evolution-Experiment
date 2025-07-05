
set -e
cd ..

VCF_DIR="results/vcf"
REF="data/reference/rel606.fna"
CHROM="NC_012967.1"
START=4375567
END=4377414
OUTDIR="results/mutLvar"
mkdir -p "$OUTDIR"

for vcf in ${VCF_DIR}/*_final_variants.vcf; do
    sample=$(basename "$vcf" _final_variants.vcf)

    echo  "Processing sample: $sample"

    # Step 1: bgzip + index if needed
    if [[ ! -f "${vcf}.gz" ]]; then
        echo "Compressing $vcf ..."
        bgzip -c "$vcf" > "${vcf}.gz"
    fi

    if [[ ! -f "${vcf}.gz.tbi" ]]; then
        echo "Indexing ${vcf}.gz ..."
        tabix -p vcf "${vcf}.gz"
    fi

    # Step 2: Extract mutL region
    mutl_out="${OUTDIR}/${sample}_mutL.vcf"
    echo "Extracting mutL variants..."
    bcftools view -r ${CHROM}:${START}-${END} "${vcf}.gz" -o "$mutl_out"



done

echo  " All mutL variant extractions complete"