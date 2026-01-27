#!/bin/bash

#SBATCH --job-name=SHAPEIT_ARRAY_TASK_WRAPPER
#SBATCH --output=SHAPEIT_ARRAY_TASK_WRAPPER.out
#SBATCH --error=SHAPEIT_ARRAY_TASK_WRAPPERsqueue.err
#SBATCH --mem-per-cpu=16000
#SBATCH --account=mignot
#SBATCH --time=120:00:00

# Parameters: $1 = prefix, $2 = dependency job ID
PREFIX=$1
DEPENDENCY=$2

# Get the pipeline directory - try multiple methods
# Method 1: Try to get it from the script location (works when script is called directly)
if [ -n "${BASH_SOURCE[0]}" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    SCRIPT_PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
else
    # Method 2: Try to get from current working directory
    SCRIPT_PIPELINE_DIR="$(pwd)"
fi

# Try to read from settings.json
SETTINGS=""
for test_dir in "$SCRIPT_PIPELINE_DIR" "$(pwd)" "/labs/mignot/researchers/smuniz/CASPR2_GWAS/imputePipeline_v2.0.0"; do
    if [ -f "${test_dir}/settings.json" ]; then
        SETTINGS="${test_dir}/settings.json"
        PIPELINE_DIR=$(jq -r '.folder.FILESFOLDER' "$SETTINGS" | sed 's|/$||')
        GWAS_BY_CHR_FOLDER=$(jq -r '.folder.GWAS_BY_CHR' "$SETTINGS")
        break
    fi
done

if [ -z "$SETTINGS" ]; then
    echo "ERROR: Could not find settings.json"
    exit 1
fi

# Change to pipeline directory to ensure we're in the right place
cd "${PIPELINE_DIR}" || { echo "ERROR: Could not cd to ${PIPELINE_DIR}"; exit 1; }

# Use GWAS_BY_CHR_FOLDER for input and output paths
INPUT_PATH="${GWAS_BY_CHR_FOLDER}"
OUTPUT_PATH="${GWAS_BY_CHR_FOLDER}"

touch shapeit_array.sh
chmod 755 shapeit_array.sh
if [[ $DEPENDENCY -eq 0 ]]; then
cat > shapeit_array.sh <<- EOF
#!/bin/bash -l
#SBATCH --job-name=shapeit_array
#SBATCH --mem-per-cpu=10000
#SBATCH --time=120:00:00
#SBATCH --array=1-22
#SBATCH --cpus-per-task=4
#SBATCH --account=mignot
cd ${PIPELINE_DIR} && ./bin/shapeit --input-bed ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.bed ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.bim ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.fam -M /labs/mignot/raw_data/gwas/1000Genomes/IMPUTE_REFERENCE_PHASE3/genetic_map_chr\${SLURM_ARRAY_TASK_ID}_combined_b37.txt -O ${OUTPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID} -T 8
EOF
else
cat > shapeit_array.sh <<- EOF
#!/bin/bash -l
#SBATCH --job-name=shapeit_array
#SBATCH --mem-per-cpu=10000
#SBATCH --time=120:00:00
#SBATCH --array=1-22
#SBATCH --depend=afterok:${DEPENDENCY}
#SBATCH --cpus-per-task=4
#SBATCH --account=mignot
cd ${PIPELINE_DIR} && ./bin/shapeit --input-bed ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.bed ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.bim ${INPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID}.fam -M /labs/mignot/raw_data/gwas/1000Genomes/IMPUTE_REFERENCE_PHASE3/genetic_map_chr\${SLURM_ARRAY_TASK_ID}_combined_b37.txt -O ${OUTPUT_PATH}${PREFIX}_CHR\${SLURM_ARRAY_TASK_ID} -T 8
EOF
fi
# Submit the array job and capture the job ID
ARRAY_JOB_ID=$(sbatch --parsable shapeit_array.sh)
echo "$ARRAY_JOB_ID"