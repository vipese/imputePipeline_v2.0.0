#!/bin/bash
###################################################################
# Batch Imputation Submission Script
# 
# This script:
# 1. Reads source path for non-imputed PLINK files from settings.json
#    (uses folder.SOURCE_DATA if set, otherwise folder.FILESFOLDER)
# 2. Reads the destination (imputed) path from settings.json (folder.BIN_FOLDER)
# 3. Identifies PLINK files that have not been imputed yet
# 4. Submits a SLURM job for each unimputed file using main.sh
#
# Each job runs in its own isolated working directory to prevent
# conflicts between concurrent jobs.
#
# Usage: ./submit_batch_imputation.sh [--dry-run] [--test]
#   --dry-run: Show what would be submitted without actually submitting
#   --test:    Only process the first 3 files (for testing)
###################################################################

set -e

# Get the directory where this script is located (same as pipeline directory)
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Parse arguments
DRY_RUN=false
TEST_MODE=false
for arg in "$@"; do
    case $arg in
        --dry-run)
            DRY_RUN=true
            echo "=== DRY RUN MODE - No jobs will be submitted ==="
            ;;
        --test)
            TEST_MODE=true
            echo "=== TEST MODE - Only processing first 3 files ==="
            ;;
    esac
done
echo ""

# Read settings from settings.json
SETTINGS="${PIPELINE_DIR}/settings.json"

if [ ! -f "$SETTINGS" ]; then
    echo "ERROR: settings.json not found at ${SETTINGS}"
    exit 1
fi

# Read paths from settings.json
# SOURCE_DATA is optional - if not set, falls back to FILESFOLDER
SOURCE_FOLDER=$(jq -r '.folder.SOURCE_DATA // .folder.FILESFOLDER' "$SETTINGS")
DEST_FOLDER=$(jq -r '.folder.BIN_FOLDER' "$SETTINGS")

# Remove trailing slash if present
SOURCE_FOLDER="${SOURCE_FOLDER%/}"
DEST_FOLDER="${DEST_FOLDER%/}"

echo "=========================================="
echo "=== Batch Imputation Submission Script ==="
echo "=========================================="
echo "Pipeline directory: $PIPELINE_DIR"
echo "Source folder (non-imputed): $SOURCE_FOLDER"
echo "Destination folder (imputed): $DEST_FOLDER"
echo ""

# Verify source folder exists
if [ ! -d "$SOURCE_FOLDER" ]; then
    echo "ERROR: Source folder does not exist: $SOURCE_FOLDER"
    exit 1
fi

# Create destination folder if it doesn't exist
if [ ! -d "$DEST_FOLDER" ]; then
    echo "Creating destination folder: $DEST_FOLDER"
    mkdir -p "$DEST_FOLDER"
fi

# Create a log directory for batch job outputs
BATCH_LOG_DIR="${PIPELINE_DIR}/BATCH_IMPUTE_LOGS"
mkdir -p "$BATCH_LOG_DIR"

# Create a base directory for job-specific working directories
WORK_BASE_DIR="${PIPELINE_DIR}/BATCH_WORKDIRS"
mkdir -p "$WORK_BASE_DIR"

# Find all PLINK files in source folder (look for .bed files as primary indicator)
echo "Scanning for PLINK files in source folder..."
echo ""

# Array to store files to impute
declare -a FILES_TO_IMPUTE

# Find all .bed files in source folder (non-recursive, just the main folder)
for bed_file in "${SOURCE_FOLDER}"/*.bed; do
    # Skip if no matches
    [ -e "$bed_file" ] || continue
    
    # Extract prefix (filename without .bed extension)
    prefix=$(basename "$bed_file" .bed)
    
    # Skip if this looks like a chromosome-split file (contains _CHR)
    if [[ "$prefix" == *"_CHR"* ]]; then
        echo "  Skipping chromosome-split file: $prefix"
        continue
    fi
    
    # Check if corresponding .bim and .fam files exist
    bim_file="${SOURCE_FOLDER}/${prefix}.bim"
    fam_file="${SOURCE_FOLDER}/${prefix}.fam"
    
    if [ ! -f "$bim_file" ]; then
        echo "  WARNING: Missing .bim file for $prefix, skipping"
        continue
    fi
    
    if [ ! -f "$fam_file" ]; then
        echo "  WARNING: Missing .fam file for $prefix, skipping"
        continue
    fi
    
    # Check if this file has already been imputed
    # Imputed files should exist as ${DEST_FOLDER}/${prefix}.bed (merged output)
    imputed_bed="${DEST_FOLDER}/${prefix}.bed"
    
    if [ -f "$imputed_bed" ]; then
        echo "  Already imputed: $prefix"
    else
        echo "  To impute: $prefix"
        FILES_TO_IMPUTE+=("$prefix")
    fi
done

echo ""
echo "=========================================="
echo "Summary:"
echo "  Total PLINK files found to impute: ${#FILES_TO_IMPUTE[@]}"

# In test mode, limit to first 3 files
if [ "$TEST_MODE" = true ] && [ ${#FILES_TO_IMPUTE[@]} -gt 3 ]; then
    FILES_TO_IMPUTE=("${FILES_TO_IMPUTE[@]:0:3}")
    echo "  Limited to first 3 files (test mode): ${#FILES_TO_IMPUTE[@]}"
fi

echo "  Files to process: ${#FILES_TO_IMPUTE[@]}"
echo "=========================================="
echo ""

if [ ${#FILES_TO_IMPUTE[@]} -eq 0 ]; then
    echo "No files need imputation. Exiting."
    exit 0
fi

# Submit jobs for each file
echo "Submitting jobs..."
echo ""

SUBMITTED_JOBS=()

for prefix in "${FILES_TO_IMPUTE[@]}"; do
    echo "Processing: $prefix"
    
    # Create a unique job name based on prefix
    JOB_NAME="IMPUTE_${prefix}"
    
    # Create unique log file paths (with timestamp to prevent overwrites)
    TIMESTAMP=$(date +%Y%m%d_%H%M%S)
    LOG_FILE="${BATCH_LOG_DIR}/${prefix}_${TIMESTAMP}.out"
    ERR_FILE="${BATCH_LOG_DIR}/${prefix}_${TIMESTAMP}.err"
    
    # Create a job-specific working directory
    # This ensures jobs don't interfere with each other
    JOB_WORKDIR="${WORK_BASE_DIR}/${prefix}_${TIMESTAMP}"
    
    if [ "$DRY_RUN" = true ]; then
        echo "  Would create workdir: $JOB_WORKDIR"
        echo "  Would submit job: $JOB_NAME"
        echo "  Output log: $LOG_FILE"
        echo "  Error log: $ERR_FILE"
        echo ""
    else
        # Create the job-specific working directory
        mkdir -p "$JOB_WORKDIR"
        
        # Copy the pipeline scripts to the working directory (symlink to save space)
        ln -sf "${PIPELINE_DIR}/scripts" "${JOB_WORKDIR}/scripts"
        ln -sf "${PIPELINE_DIR}/bin" "${JOB_WORKDIR}/bin"
        cp "${PIPELINE_DIR}/main.sh" "${JOB_WORKDIR}/main.sh"
        
        # Create job-specific settings.json with the correct prefix and paths
        # Use the job workdir as FILESFOLDER since data will be copied there
        jq --arg prefix "$prefix" \
           --arg filesfolder "${JOB_WORKDIR}/" \
           --arg gwas_by_chr "${JOB_WORKDIR}/GWAS_BY_CHR/" \
           --arg slurm_log "${JOB_WORKDIR}/SLURM_IMPUTE_LOG/" \
           --arg shapeit_log "${JOB_WORKDIR}/SHAPEIT_IMPUTE_LOG/" \
           '.prefix = $prefix | 
            .folder.FILESFOLDER = $filesfolder |
            .folder.GWAS_BY_CHR = $gwas_by_chr |
            .folder.SLURM_IMPUTE_LOG = $slurm_log |
            .folder.SHAPEIT_IMPUTE_LOG = $shapeit_log' \
           "$SETTINGS" > "${JOB_WORKDIR}/settings.json"
        
        # Copy source PLINK files to the job working directory
        cp "${SOURCE_FOLDER}/${prefix}.bed" "${JOB_WORKDIR}/"
        cp "${SOURCE_FOLDER}/${prefix}.bim" "${JOB_WORKDIR}/"
        cp "${SOURCE_FOLDER}/${prefix}.fam" "${JOB_WORKDIR}/"
        
        # Submit the job
        JOB_ID=$(sbatch --parsable \
            --job-name="$JOB_NAME" \
            --output="$LOG_FILE" \
            --error="$ERR_FILE" \
            --mem-per-cpu=16000 \
            --time=48:00:00 \
            --account=mignot \
            --wrap="
                set -e
                
                # Change to job-specific working directory
                cd ${JOB_WORKDIR}
                
                echo '=========================================='
                echo 'Job: ${JOB_NAME}'
                echo 'Working directory: ${JOB_WORKDIR}'
                echo 'Prefix: ${prefix}'
                echo 'Started at: \$(date)'
                echo '=========================================='
                echo ''
                
                # Run the main imputation pipeline
                bash main.sh
                
                echo ''
                echo '=========================================='
                echo 'Moving results to destination folder...'
                echo '=========================================='
                
                # Move the final imputed results to the destination folder
                # The merged output should be in BIN_FOLDER as specified in settings
                BIN_FOLDER=\$(jq -r '.folder.BIN_FOLDER' settings.json)
                
                if [ -d \"\$BIN_FOLDER\" ]; then
                    # Move all output files to the final destination
                    cp -r \"\$BIN_FOLDER\"/* \"${DEST_FOLDER}/\" 2>/dev/null || true
                    echo 'Results copied to ${DEST_FOLDER}'
                fi
                
                echo ''
                echo '=========================================='
                echo 'Cleaning up working directory...'
                echo '=========================================='
                
                # Clean up the job working directory to save space
                cd ${WORK_BASE_DIR}
                rm -rf ${JOB_WORKDIR}
                
                echo 'Imputation completed for ${prefix} at \$(date)'
            ")
        
        if [ -n "$JOB_ID" ]; then
            echo "  Created workdir: $JOB_WORKDIR"
            echo "  Submitted job $JOB_ID: $JOB_NAME"
            SUBMITTED_JOBS+=("$JOB_ID")
        else
            echo "  ERROR: Failed to submit job for $prefix"
            # Clean up the working directory if job submission failed
            rm -rf "$JOB_WORKDIR"
        fi
        echo ""
    fi
done

echo "=========================================="
echo "=== Submission Complete ==="
echo "=========================================="

if [ "$DRY_RUN" = false ]; then
    echo "Submitted ${#SUBMITTED_JOBS[@]} jobs"
    echo ""
    echo "Job IDs: ${SUBMITTED_JOBS[*]}"
    echo ""
    echo "Log files are in: $BATCH_LOG_DIR"
    echo "Working directories are in: $WORK_BASE_DIR"
    echo ""
    echo "To monitor jobs, run:"
    echo "  squeue -u \$USER"
    echo ""
    echo "To check a specific job:"
    echo "  squeue -j <job_id>"
    echo ""
    echo "Note: Working directories will be automatically cleaned up after each job completes."
else
    echo ""
    echo "Run without --dry-run to actually submit jobs."
fi
