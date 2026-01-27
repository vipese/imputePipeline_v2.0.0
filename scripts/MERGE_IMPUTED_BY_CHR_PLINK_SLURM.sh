#!/bin/bash -l
#SBATCH --job-name=MERGE_BY_CHR
#SBATCH --mem-per-cpu=64000
#SBATCH --account=mignot
#SBATCH --time=12:00:00

# Load module 
module load plink

# Go to path 
path=/labs/mignot/Lyon_LGI1_Imputed/
cd $path 

# Create mergelist 
find ./*.bed -printf "%f\n" | cut -d '.' -f 1 > mergeList.txt

# Initialize plink input
CHRFILES=($(awk '{print $0}' mergeList.txt))
LENFILES=${#CHRFILES[@]}

# Perform QC
for IDX in $(seq 0 $LENFILES); do

   # Parse duplicate variants (based on position and allele codes)
    plink --bfile ${CHRFILES[IDX]} --list-duplicate-vars suppress-first \
        --out temp > temp
    awk '{print $4}' temp.dupvar > ${CHRFILES[IDX]}.dupvar
    rm -r temp*

    echo "---"

    # Perform Quality control - Remove duplicated variants
    plink --bfile ${CHRFILES[IDX]} --exclude ${CHRFILES[IDX]}.dupvar \
        --no-sex --no-parents --not-chr 25,26 \
        --maf 0.05 \
        --make-bed --out ${CHRFILES[IDX]}_tmpFile

    echo "-----------------------------------------------------------------"
    echo "-----------------------------------------------------------------"

done

# Modify mergelist to create plink version of mergelist 
awk '{print $0 "_tmpFile.bed" " " $0 "_tmpFile.bim" " " $0 "_tmpFile.fam"}' mergeList.txt > temp && mv temp mergeList.txt

# Merge
INITFILE=$(awk 'FNR==1 {print $1}' mergeList.txt | cut -d '.' -f 1)
plink --bfile $INITFILE --merge-list mergeList.txt --out Lyon_LGI1

# Remove temporal files 
rm -r *tmpFile*