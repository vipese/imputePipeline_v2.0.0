#!/bin/sh
#1 - chromosome number
#2 total length of chr
#3 file name
#4 shapeit_jobID

# Get the absolute path of the pipeline directory (where this script lives)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
cd "$PIPELINE_DIR"

mkdir -p imputeFiles

# Get the first SNP position from the .haps file and convert to MB (rounded down)
# Check in current directory first, then in GWAS_BY_CHR/
HAPS_FILE="$3"_CHR"$1".haps
if [ ! -f "$HAPS_FILE" ]; then
    # Try GWAS_BY_CHR directory
    if [ -f "GWAS_BY_CHR/$HAPS_FILE" ]; then
        HAPS_FILE="GWAS_BY_CHR/$HAPS_FILE"
    else
        echo "Error: $HAPS_FILE not found in current directory or GWAS_BY_CHR/!"
        echo "Current directory: $(pwd)"
        echo "Looking for: $HAPS_FILE"
        exit 1
    fi
fi

# Extract first SNP position (column 3), divide by 1000000 to get MB, round down
START_MB=$(awk 'NR==1 {print int($3/1000000)}' "$HAPS_FILE")
echo "Chromosome $1: Starting imputation from ${START_MB}MB (first SNP detected)"

# Start loop from first SNP position instead of 0
for i in `seq $START_MB $2` 
do
interval="${i}e6 $(($i + 1))e6"

# Use absolute path for the working directory
WORKDIR="$PIPELINE_DIR"
LOGDIR="${WORKDIR}/SLURM_IMPUTE_LOG"
mkdir -p "$LOGDIR"

# Build command with absolute paths
IMPUTE_CMD="${WORKDIR}/bin/impute2 -known_haps_g ${WORKDIR}/${HAPS_FILE} -h /labs/mignot/raw_data/gwas/1000Genomes/IMPUTE_REFERENCE_PHASE3/1000GP_Phase3_chr${1}.hap.gz -l /labs/mignot/raw_data/gwas/1000Genomes/IMPUTE_REFERENCE_PHASE3/1000GP_Phase3_chr${1}.legend.gz -m /labs/mignot/raw_data/gwas/1000Genomes/IMPUTE_REFERENCE_PHASE3/genetic_map_chr${1}_combined_b37.txt -int ${interval} -buffer 500 -Ne 20000 -o ${WORKDIR}/imputeFiles/CHR${1}_${3}.${i}"

if [[ $4 -eq 0 ]]; then
cat > tmpchr"$1".$i.sh <<-EOF
#!/bin/bash -l
#SBATCH --job-name=EM_$i.chr"$1"
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH --account=mignot
#SBATCH --output=${LOGDIR}/impute_chr${1}_${i}_%j.out
#SBATCH --error=${LOGDIR}/impute_chr${1}_${i}_%j.err
$IMPUTE_CMD 
EOF
else
cat > tmpchr"$1".$i.sh <<-EOF
#!/bin/bash -l
#SBATCH --job-name=EM_$i.chr"$1"
#SBATCH --depend=afterok:"$4"_"$1"
#SBATCH --mem-per-cpu=10000
#SBATCH --time=12:00:00
#SBATCH --account=mignot
#SBATCH --output=${LOGDIR}/impute_chr${1}_${i}_%j.out
#SBATCH --error=${LOGDIR}/impute_chr${1}_${i}_%j.err
$IMPUTE_CMD
EOF
fi
pending=$(squeue -t pd -u $USER -h | wc -l)
sbatch tmpchr"$1".$i.sh
while [[ ${pending} -gt 100 ]]
do
sleep 60
pending=$(squeue -t pd -u $USER -h | wc -l)
done
rm tmpchr"$1".$i.sh
done

