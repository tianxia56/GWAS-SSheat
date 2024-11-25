#!/bin/bash

#!/bin/bash

# Function to run plink clumping
run_plink_clump() {
    local bfile=$1
    local input_file=$2
    local output_prefix=$3

    plink --bfile "$bfile" --clump-p1 1 --clump-r2 0.1 --clump-kb 250 \
          --clump "$input_file" --clump-snp-field SNP --clump-field P \
          --out "$output_prefix"
}

# Define the combinations
declare -A files
files=( ["Height.BBJ.QC_with_SNP.txt"]="jpt.all" 
        ["Height.UKBB.QC_with_SNP.txt"]="gbr.all" 
        ["TG.BBJ.QC_with_SNP.txt"]="jpt.all" 
        ["TG.UKBB.QC_with_SNP.txt"]="gbr.all" 
        ["BMI.BBJ.QC_with_SNP.txt"]="jpt.all" 
        ["BMI.UKBB.QC_with_SNP.txt"]="gbr.all" )

# Run plink clumping for each combination
for file in "${!files[@]}"; do
    bfile=${files[$file]}
    output_prefix=$(basename "$file" .QC_with_SNP.txt)
    run_plink_clump "$bfile" "$file" "$output_prefix"
done
