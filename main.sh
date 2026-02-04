#!/bin/bash -l
#SBATCH --job-name=FULL_IMPUTE_PIPELINE
#SBATCH --mem-per-cpu=16000
#SBATCH --time=24:00:00
#SBATCH --account=mignot
#SBATCH --output=FULL_IMPUTE_PIPELINE.out
#SBATCH --error=FULL_IMPUTE_PIPELINE.err

###################################################################
# Full End-to-End Imputation Pipeline
# Combines preprocessing, imputation, and post-imputation steps
# Based on main.sh with fixes from run_imputation.sh and run_post_imputation.sh
###################################################################

# cd /labs/mignot/researchers/smuniz/CASPR2_GWAS/imputePipeline_v2.0.0

# Read settings
SETTINGS=$(pwd)/settings.json
PREFIX=$(jq -r '.prefix' $SETTINGS)
REF=$(jq -r '.ref' $SETTINGS)

# Cleanup option: set to true to remove intermediate files at the end
CLEANUP=true

# Optional: set to true to remove per-chromosome files from destination folder after merge
# This removes CHR*_${PREFIX}.* files (BGEN, BED, etc.) leaving only the merged output
CLEANUP_CHR_FILES=false

# Initialize folders
FILESFOLDER=$(jq -r '.folder.FILESFOLDER' $SETTINGS)
GWAS_BY_CHR_FOLDER=$(jq -r '.folder.GWAS_BY_CHR' $SETTINGS)
SLURM_IMPUTE_LOG=$(jq -r '.folder.SLURM_IMPUTE_LOG' $SETTINGS)
SHAPEIT_IMPUTE_LOG=$(jq -r '.folder.SHAPEIT_IMPUTE_LOG' $SETTINGS)
BINFILES_FOLDER=$(jq -r '.folder.BIN_FOLDER' $SETTINGS)
SCRIPTS=${FILESFOLDER}scripts/

echo "=========================================="
echo "=== Full Imputation Pipeline Starting ==="
echo "=========================================="
echo "PREFIX: $PREFIX"
echo "Reference: $REF"
echo "Date: $(date)"
echo ""

########## STEP 1: PREPROCESSING ###########
echo "=== Step 1: Data Preprocessing ==="

# Load module
module load plink

# If files are not binarized, binarize
if [ ! -f ${PREFIX}.bed ]; then
    echo "Binarizing PLINK files..."
    plink --file ${PREFIX} --allow-no-sex --no-sex --no-fid --no-parents --no-pheno --make-bed --out ${PREFIX}
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to binarize PLINK files"
        exit 1
    fi
else
    echo "PLINK binary files already exist. Skipping binarization."
fi

# Create list of duplicated variants by name and position
echo "Identifying duplicate variants..."
awk 'a[$2]++{print $0}' ${PREFIX}.bim > Dup_vars_name_pos.txt
awk 'a[$4]++{print $0}' ${PREFIX}.bim >> Dup_vars_name_pos.txt
sort -k3n ${PREFIX}.bim | uniq -f2 -D >> Dup_vars_name_pos.txt

# Iteratively remove multi-allelic variants (PLINK does not manage them well)
NDUPVARS=$(awk 'END{print NR}' Dup_vars_name_pos.txt)
ITER=0
while [ $NDUPVARS -gt 1 ]; do
    ITER=$((ITER + 1))
    echo "Removing duplicates - iteration $ITER (found $NDUPVARS duplicates)..."
    
    plink --bfile ${PREFIX} --exclude Dup_vars_name_pos.txt \
        --allow-no-sex \
        --make-bed --out gwastempFilt > gwastempFilt.log 2>&1
    
    if [ $? -ne 0 ]; then
        echo "ERROR: Failed to remove duplicates in iteration $ITER"
        exit 1
    fi
    
    rm -f ${PREFIX}.bed ${PREFIX}.fam ${PREFIX}.bim
    mv gwastempFilt.bed ${PREFIX}.bed
    mv gwastempFilt.fam ${PREFIX}.fam
    mv gwastempFilt.bim ${PREFIX}.bim
    rm -f gwastempFilt*

    awk 'a[$2]++{print $0}' ${PREFIX}.bim > Dup_vars_name_pos.txt
    sort -k3n ${PREFIX}.bim | uniq -f2 -D >> Dup_vars_name_pos.txt
    NDUPVARS=$(awk 'END{print NR}' Dup_vars_name_pos.txt)
done
echo "Duplicate removal completed."

# Filter out low genotyping in samples and variants (SHAPEIT throws an error at fully missing vars or subjects)
echo "Filtering low genotyping samples and variants..."
plink --bfile ${PREFIX} --allow-no-sex \
    --mind 0.05 \
    --geno 0.05 \
    --make-bed --out gwastempFilt > gwastempFilt.log 2>&1

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to filter low genotyping"
    exit 1
fi

mv gwastempFilt.bed ${PREFIX}.bed
mv gwastempFilt.fam ${PREFIX}.fam
mv gwastempFilt.bim ${PREFIX}.bim
rm -f gwastempFilt*
echo "Preprocessing completed at $(date)"
echo ""

########## STEP 2: SPLIT BY CHROMOSOME ##########
echo "=== Step 2: Splitting by Chromosome ==="

# Check if chromosome-split files exist
CHR_SPLIT_EXISTS=true
for chr in {1..22}; do
    if [ ! -f "${PREFIX}_CHR${chr}.bed" ] && [ ! -f "${GWAS_BY_CHR_FOLDER}${PREFIX}_CHR${chr}.bed" ]; then
        CHR_SPLIT_EXISTS=false
        break
    fi
done

if [ "$CHR_SPLIT_EXISTS" = false ]; then
    echo "Chromosome-split files not found. Running PLINK split..."
    
    # Submit PLINK split job
    SPLIT_JOB=$(sbatch --parsable --time=01:00:00 ${SCRIPTS}PLINK_SPLIT_SLURM.sh $PREFIX)
    if [ -z "$SPLIT_JOB" ]; then
        echo "ERROR: Failed to submit PLINK split job"
        exit 1
    fi
    echo "PLINK split job: $SPLIT_JOB"
    
    # Wait for split to complete
    sleep 60  # Initial wait for jobs to appear
    while squeue -j $SPLIT_JOB -h 2>/dev/null | grep -q .; do
        RUNNING=$(squeue -j $SPLIT_JOB -h 2>/dev/null | grep -c "R" || echo 0)
        PENDING=$(squeue -j $SPLIT_JOB -h 2>/dev/null | grep -c "PD" || echo 0)
        echo "$(date): PLINK split - Running: $RUNNING, Pending: $PENDING"
        sleep 300
    done
    echo "PLINK split completed at $(date)"
    
    # Move split files to GWAS_BY_CHR if they're in current directory
    mkdir -p ${GWAS_BY_CHR_FOLDER}
    for chr in {1..22}; do
        if [ -f "${PREFIX}_CHR${chr}.bed" ]; then
            mv ${PREFIX}_CHR${chr}.* ${GWAS_BY_CHR_FOLDER}/ 2>/dev/null || true
        fi
    done
    
    # Verify split files exist
    SPLIT_COUNT=0
    for chr in {1..22}; do
        if [ -f "${GWAS_BY_CHR_FOLDER}${PREFIX}_CHR${chr}.bed" ]; then
            SPLIT_COUNT=$((SPLIT_COUNT + 1))
        fi
    done
    
    if [ "$SPLIT_COUNT" -lt 22 ]; then
        echo "ERROR: Only found $SPLIT_COUNT chromosome-split files. Expected 22."
        exit 1
    fi
    echo "All 22 chromosome-split files created."
else
    echo "Chromosome-split files already exist. Skipping split step."
fi
echo ""

########## STEP 3: PHASING WITH SHAPEIT2 ##########
echo "=== Step 3: Phasing with SHAPEIT2 ==="

# Check if .haps files exist
HAPS_EXISTS=false
for chr in {1..22}; do
    if [ ! -f "${PREFIX}_CHR${chr}.haps" ] && [ ! -f "${GWAS_BY_CHR_FOLDER}${PREFIX}_CHR${chr}.haps" ]; then
        HAPS_EXISTS=false
        break
    fi
done

if [ "$HAPS_EXISTS" = false ]; then
    echo ".haps files not found. Running SHAPEIT2 phasing..."
    
    # Check if we have a dependency from split job
    if [ -n "$SPLIT_JOB" ]; then
        SHAPEIT_JOB=$(sbatch --parsable --time=120:00:00 ${SCRIPTS}SHAPEIT_ARRAY_TASK_SLURM.sh $PREFIX $SPLIT_JOB)
    else
        SHAPEIT_JOB=$(sbatch --parsable --time=120:00:00 ${SCRIPTS}SHAPEIT_ARRAY_TASK_SLURM.sh $PREFIX 0)
    fi
    
    if [ -z "$SHAPEIT_JOB" ]; then
        echo "ERROR: Failed to submit SHAPEIT2 phasing job"
        exit 1
    fi
    echo "SHAPEIT2 phasing job: $SHAPEIT_JOB"
    
    # Wait for phasing to complete
    sleep 60  # Initial wait for jobs to appear
    while true; do
        RUNNING=$(squeue -u $USER -h -n "shapeit_array" 2>/dev/null | wc -l)
        PENDING=$(squeue -u $USER -h -t PD 2>/dev/null | grep "shapeit_array" | wc -l)
        
        if [ "$RUNNING" -eq 0 ] && [ "$PENDING" -eq 0 ]; then
            sleep 30
            RUNNING=$(squeue -u $USER -h 2>/dev/null | grep "shapeit_array" | wc -l)
            if [ "$RUNNING" -eq 0 ]; then
                echo "SHAPEIT2 phasing completed at $(date)"
                break
            fi
        fi
        
        echo "$(date): SHAPEIT2 - Running: $RUNNING, Pending: $PENDING"
        sleep 300
    done
    
    # Move .haps files to GWAS_BY_CHR if needed
    mkdir -p ${GWAS_BY_CHR_FOLDER}
    for chr in {1..22}; do
        if [ -f "${PREFIX}_CHR${chr}.haps" ]; then
            mv ${PREFIX}_CHR${chr}.haps ${PREFIX}_CHR${chr}.haps.sample ${GWAS_BY_CHR_FOLDER}/ 2>/dev/null || true
        fi
    done
    
    # Verify .haps files exist
    HAPS_COUNT=0
    for chr in {1..22}; do
        if [ -f "${PREFIX}_CHR${chr}.haps" ] || [ -f "${GWAS_BY_CHR_FOLDER}${PREFIX}_CHR${chr}.haps" ]; then
            HAPS_COUNT=$((HAPS_COUNT + 1))
        fi
    done
    
    if [ "$HAPS_COUNT" -lt 22 ]; then
        echo "ERROR: Only found $HAPS_COUNT .haps files. Expected 22. Check SHAPEIT logs."
        exit 1
    fi
    echo "All 22 .haps files created."
else
    echo ".haps files already exist. Skipping phasing step."
fi
echo ""

########## STEP 4: IMPUTATION WITH IMPUTE2 ##########
echo "=== Step 4: Imputation with IMPUTE2 ==="

# Chromosome sizes for imputation (in MB)
declare -A chrSizes=(
    [1]=250 [2]=244 [3]=199 [4]=192 [5]=181 [6]=172 [7]=160 [8]=147 [9]=142
    [10]=136 [11]=136 [12]=134 [13]=116 [14]=108 [15]=103 [16]=91 [17]=82
    [18]=79 [19]=60 [20]=64 [21]=49 [22]=52
)

# Verify .haps files exist before submitting imputation
echo "Verifying .haps files exist..."
MISSING_HAPS=0
for chr in {1..22}; do
    if [ ! -f "${PREFIX}_CHR${chr}.haps" ] && [ ! -f "${GWAS_BY_CHR_FOLDER}${PREFIX}_CHR${chr}.haps" ]; then
        echo "ERROR: ${PREFIX}_CHR${chr}.haps not found!"
        MISSING_HAPS=$((MISSING_HAPS + 1))
    fi
done

if [ "$MISSING_HAPS" -gt 0 ]; then
    echo "ERROR: Missing $MISSING_HAPS .haps files. Cannot proceed with imputation."
    exit 1
fi

# Create imputeFiles directory
mkdir -p imputeFiles

# Submit all imputation jobs
echo "Submitting IMPUTE2 jobs for all chromosomes..."
for chr in {1..22}; do
    chrSize=${chrSizes[$chr]}
    echo "Submitting CHR$chr (size: ${chrSize}MB)..."
    
    # Run the impute loop script - it will submit individual segment jobs
    if ! ./scripts/IMPUTE_LOOP_SLURM_FIXED.sh $chr $chrSize $PREFIX 0; then
        echo "ERROR: Failed to submit imputation jobs for CHR$chr"
        exit 1
    fi
done

echo "All imputation jobs submitted."
echo "Waiting for imputation jobs to complete (this may take many hours)..."

# Wait for all EM_ jobs to complete
sleep 60  # Initial wait for jobs to appear in queue
while true; do
    RUNNING=$(squeue -u $USER -h -n "EM_*" 2>/dev/null | wc -l)
    PENDING=$(squeue -u $USER -h -t PD 2>/dev/null | grep "EM_" | wc -l)
    
    if [ "$RUNNING" -eq 0 ] && [ "$PENDING" -eq 0 ]; then
        # Double check after a short wait
        sleep 30
        RUNNING=$(squeue -u $USER -h 2>/dev/null | grep "EM_" | wc -l)
        if [ "$RUNNING" -eq 0 ]; then
            echo "All imputation jobs completed at $(date)"
            break
        fi
    fi
    
    echo "$(date): IMPUTE2 - Running: $RUNNING, Pending: $PENDING"
    sleep 300  # Check every 5 minutes
done

# Verify imputation output exists
echo "Checking for imputation output files..."
IMPUTE_COUNT=$(find imputeFiles/ -name "CHR*_${PREFIX}.*" -type f 2>/dev/null | wc -l)
echo "Found $IMPUTE_COUNT total imputation output files"

# Check for actual impute output files (not just info files)
IMPUTE_OUTPUT_COUNT=$(find imputeFiles/ -name "CHR*_${PREFIX}.[0-9]*" -type f 2>/dev/null | wc -l)
echo "Found $IMPUTE_OUTPUT_COUNT imputation segment files"

if [ "$IMPUTE_OUTPUT_COUNT" -lt 100 ]; then
    echo "ERROR: Expected at least 100 imputation segment files, but found only $IMPUTE_OUTPUT_COUNT"
    echo "Check logs in SLURM_IMPUTE_LOG/ for errors"
    
    # Check for common error patterns in logs
    if [ -d "SLURM_IMPUTE_LOG" ]; then
        echo "Recent errors from imputation logs:"
        find SLURM_IMPUTE_LOG/ -name "*.err" -type f -exec tail -5 {} \; 2>/dev/null | head -20
    fi
    exit 1
fi
echo ""

########## STEP 5: CLEAN UP INTERMEDIATE FILES ##########
echo "=== Step 5: Cleaning up intermediate files ==="

# Clean up: move CHR files to directory (but keep .fam files in original location)
# Note: .fam files are NOT moved because:
#   1. They're still needed in Step 7 (sort/convert to BGEN) to create sample files
#   2. All chromosomes have the same samples, so we can use the original ${PREFIX}.fam
#   3. This avoids path issues and keeps files where they're expected
mkdir -p $GWAS_BY_CHR_FOLDER
# Move all CHR files EXCEPT .fam files (they stay in original location)
for file in ${PREFIX}_CHR*.*; do
    if [ -f "$file" ] && [[ "$file" != *.fam ]]; then
        mv "$file" $GWAS_BY_CHR_FOLDER/ 2>/dev/null || true
    fi
done

# Clean up: move slurm outputs
mkdir -p $SLURM_IMPUTE_LOG
mv slurm-*.out $SLURM_IMPUTE_LOG/ 2>/dev/null || true
mv slurm-*.err $SLURM_IMPUTE_LOG/ 2>/dev/null || true

# Clean up shapeit logs
mkdir -p $SHAPEIT_IMPUTE_LOG
mv shapeit*.log $SHAPEIT_IMPUTE_LOG/ 2>/dev/null || true
mv shapeit*.out $SHAPEIT_IMPUTE_LOG/ 2>/dev/null || true
echo "Cleanup completed."
echo ""

########## STEP 6: CONCATENATE IMPUTED SEGMENTS ##########
echo "=== Step 6: Concatenating imputed segments ==="

# Clean BINFILES folder / create directory
if [ -d $BINFILES_FOLDER ]; then
    rm -rf $BINFILES_FOLDER
fi
mkdir -p $BINFILES_FOLDER

# Run concatenation
CAT_JOB=$(sbatch --parsable --time=12:00:00 ${SCRIPTS}CAT_IMPUTE_SLURM.sh -d $(pwd)/imputeFiles/ -s $BINFILES_FOLDER -c $SCRIPTS -p $PREFIX)
if [ -z "$CAT_JOB" ]; then
    echo "ERROR: Failed to submit concatenation job"
    exit 1
fi
echo "Concatenation job: $CAT_JOB"

# Wait for concatenation
sleep 60
while squeue -j $CAT_JOB -h 2>/dev/null | grep -q .; do
    RUNNING=$(squeue -j $CAT_JOB -h 2>/dev/null | grep -c "R" || echo 0)
    PENDING=$(squeue -j $CAT_JOB -h 2>/dev/null | grep -c "PD" || echo 0)
    echo "$(date): Concatenation - Running: $RUNNING, Pending: $PENDING"
    sleep 300
done
echo "Concatenation completed at $(date)"

# Verify concatenated files exist
CONCAT_COUNT=$(ls ${BINFILES_FOLDER}CHR*_${PREFIX}.impute.gz 2>/dev/null | wc -l)
if [ "$CONCAT_COUNT" -lt 22 ]; then
    echo "ERROR: Only found $CONCAT_COUNT concatenated files. Expected 22."
    exit 1
fi
echo "All 22 concatenated files created."
echo ""

########## STEP 7: SORT AND CONVERT TO BGEN ##########
echo "=== Step 7: Sorting and converting to BGEN ==="

SORT_JOB=$(sbatch --parsable --time=12:00:00 ${SCRIPTS}SORT_IMPUTED_SLURM.sh -d $BINFILES_FOLDER -p $PREFIX)
if [ -z "$SORT_JOB" ]; then
    echo "ERROR: Failed to submit sort job"
    exit 1
fi
echo "Sort job: $SORT_JOB"

# Wait for sort
sleep 60
while squeue -j $SORT_JOB -h 2>/dev/null | grep -q .; do
    RUNNING=$(squeue -j $SORT_JOB -h 2>/dev/null | grep -c "R" || echo 0)
    PENDING=$(squeue -j $SORT_JOB -h 2>/dev/null | grep -c "PD" || echo 0)
    echo "$(date): Sort/Convert - Running: $RUNNING, Pending: $PENDING"
    sleep 600  # Check every 10 minutes (this step takes a long time)
done
echo "Sorting and BGEN conversion completed at $(date)"

# Verify BGEN files exist and are large enough
echo "Verifying BGEN files..."
BGEN_COUNT=0
for chr in {1..22}; do
    BGEN_FILE="${BINFILES_FOLDER}CHR${chr}_${PREFIX}.bgen"
    if [ -f "$BGEN_FILE" ]; then
        SIZE=$(stat -c%s "$BGEN_FILE" 2>/dev/null || echo 0)
        if [ "$SIZE" -lt 10000000 ]; then  # Less than 10MB is suspicious
            echo "WARNING: CHR${chr} BGEN file is too small ($SIZE bytes)"
        else
            echo "CHR${chr}: $(ls -lh $BGEN_FILE | awk '{print $5}')"
            BGEN_COUNT=$((BGEN_COUNT + 1))
        fi
    else
        echo "ERROR: CHR${chr} BGEN file missing!"
    fi
done

if [ "$BGEN_COUNT" -lt 22 ]; then
    echo "ERROR: Only found $BGEN_COUNT BGEN files. Expected 22."
    exit 1
fi
echo "All 22 BGEN files created and validated."
echo ""

########## STEP 8: CONVERT TO PLINK ##########
echo "=== Step 8: Converting BGEN to PLINK ==="

BGEN_JOB=$(sbatch --parsable --time=8:00:00 ${SCRIPTS}bgen2plink.sh)
if [ -z "$BGEN_JOB" ]; then
    echo "ERROR: Failed to submit BGEN2PLINK job"
    exit 1
fi
echo "BGEN2PLINK job: $BGEN_JOB"

# Wait for conversion
sleep 60
while squeue -j $BGEN_JOB -h 2>/dev/null | grep -q .; do
    RUNNING=$(squeue -j $BGEN_JOB -h 2>/dev/null | grep -c "R" || echo 0)
    PENDING=$(squeue -j $BGEN_JOB -h 2>/dev/null | grep -c "PD" || echo 0)
    echo "$(date): BGEN2PLINK - Running: $RUNNING, Pending: $PENDING"
    sleep 300
done
echo "BGEN to PLINK conversion completed at $(date)"

# Verify BED files
echo "Verifying BED files..."
BED_COUNT=$(ls ${BINFILES_FOLDER}CHR*_${PREFIX}.bed 2>/dev/null | wc -l)
echo "Found $BED_COUNT BED files"

if [ "$BED_COUNT" -lt 22 ]; then
    echo "ERROR: Only found $BED_COUNT BED files. Expected 22."
    exit 1
fi
echo ""

########## STEP 9: MERGE CHROMOSOMES ##########
echo "=== Step 9: Merging chromosomes ==="

MERGE_JOB=$(sbatch --parsable --time=8:00:00 ${SCRIPTS}merge_cohort.sh)
if [ -z "$MERGE_JOB" ]; then
    echo "ERROR: Failed to submit merge job"
    exit 1
fi
echo "Merge job: $MERGE_JOB"

# Wait for merge
sleep 60
while squeue -j $MERGE_JOB -h 2>/dev/null | grep -q .; do
    RUNNING=$(squeue -j $MERGE_JOB -h 2>/dev/null | grep -c "R" || echo 0)
    PENDING=$(squeue -j $MERGE_JOB -h 2>/dev/null | grep -c "PD" || echo 0)
    echo "$(date): Merge - Running: $RUNNING, Pending: $PENDING"
    sleep 300
done
echo "Merge completed at $(date)"

# Verify merged files
if [ ! -f "${BINFILES_FOLDER}${PREFIX}.bed" ]; then
    echo "ERROR: Merged BED file not found!"
    exit 1
fi
echo ""

########## STEP 10: FINAL CLEANUP ##########
if [ "$CLEANUP" = true ]; then
    echo "=== Step 10: Final Cleanup ==="

    # Remove SLURM imputation logs
    if [ -d "$SLURM_IMPUTE_LOG" ]; then
        echo "Removing SLURM imputation logs: $SLURM_IMPUTE_LOG"
        rm -rf "$SLURM_IMPUTE_LOG"
    fi

    # Remove SHAPEIT logs
    if [ -d "$SHAPEIT_IMPUTE_LOG" ]; then
        echo "Removing SHAPEIT logs: $SHAPEIT_IMPUTE_LOG"
        rm -rf "$SHAPEIT_IMPUTE_LOG"
    fi

    # Remove GWAS_BY_CHR folder (intermediate chromosome-split files)
    if [ -d "$GWAS_BY_CHR_FOLDER" ]; then
        echo "Removing GWAS_BY_CHR folder: $GWAS_BY_CHR_FOLDER"
        rm -rf "$GWAS_BY_CHR_FOLDER"
    fi

    # Remove imputeFiles folder (intermediate imputation segments)
    if [ -d "imputeFiles" ]; then
        echo "Removing imputeFiles folder"
        rm -rf imputeFiles
    fi

    # Remove all .out and .err files in the pipeline directory
    echo "Removing .out and .err files in pipeline directory..."
    find . -maxdepth 1 -name "*.out" -type f -delete 2>/dev/null
    find . -maxdepth 1 -name "*.err" -type f -delete 2>/dev/null

    # Remove any remaining slurm output files
    rm -f slurm-*.out slurm-*.err 2>/dev/null

    # Remove SLURM_IMPUTE_LOG folder if it was recreated
    rm -rf SLURM_IMPUTE_LOG 2>/dev/null

    # Remove any other log directories that might have been created
    rm -rf shapeit_logs 2>/dev/null
    rm -rf impute_logs 2>/dev/null

    # Remove all shapeit-related files
    echo "Removing shapeit log files..."
    rm -f shapeit*.log shapeit*.out shapeit*.err 2>/dev/null
    find . -maxdepth 1 -name "shapeit*" -type f -delete 2>/dev/null

    # Clean up temporary files
    rm -f Dup_vars_name_pos.txt 2>/dev/null
    rm -f gwastempFilt.* 2>/dev/null
    rm -f *.log 2>/dev/null

    echo "Final cleanup completed."
    echo ""
else
    echo "=== Step 10: Final Cleanup (SKIPPED) ==="
    echo "CLEANUP is set to false. Intermediate files retained."
    echo ""
fi

# Optional: Remove per-chromosome files from destination folder
if [ "$CLEANUP_CHR_FILES" = true ]; then
    echo "=== Removing per-chromosome files from destination folder ==="
    
    # Count files before removal
    CHR_FILE_COUNT=$(ls ${BINFILES_FOLDER}CHR*_${PREFIX}.* 2>/dev/null | wc -l)
    
    if [ "$CHR_FILE_COUNT" -gt 0 ]; then
        echo "Removing $CHR_FILE_COUNT per-chromosome files..."
        rm -f ${BINFILES_FOLDER}CHR*_${PREFIX}.* 2>/dev/null
        echo "Per-chromosome files removed from: $BINFILES_FOLDER"
    else
        echo "No per-chromosome files found to remove."
    fi
    echo ""
fi

########## DONE ##########
echo "=========================================="
echo "=== Pipeline Complete ==="
echo "=========================================="
echo "Finished at $(date)"
echo ""
echo "Output files in: $BINFILES_FOLDER"
echo ""
echo "Merged files:"
ls -lh ${BINFILES_FOLDER}${PREFIX}.bed ${BINFILES_FOLDER}${PREFIX}.bim ${BINFILES_FOLDER}${PREFIX}.fam 2>/dev/null || echo "WARNING: Merged files not found"
echo ""

if [ "$CLEANUP_CHR_FILES" = true ]; then
    echo "Per-chromosome files: REMOVED (CLEANUP_CHR_FILES=true)"
else
    echo "Per-chromosome BGEN files:"
    ls -lh ${BINFILES_FOLDER}CHR*_${PREFIX}.bgen 2>/dev/null | head -5
fi
echo ""
echo "All intermediate files and logs have been cleaned up."
echo ""

