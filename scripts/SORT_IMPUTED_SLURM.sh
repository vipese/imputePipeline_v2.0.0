#!/bin/bash -l
#SBATCH --job-name=SORT_IMPUTED_SLURM``
#SBATCH --mem-per-cpu=32000
#SBATCH --output=SORT_IMPUTED_SLURM.out
#SBATCH --error=SORT_IMPUTED_SLURM.err
#SBATCH --array=1-22
#SBATCH --account=mignot
#SBATCH --time=12:00:00

# Parse arguments
PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-p <path>]

-d <directory>: Path where concatenated outputs from imputePipe.py are located
-p <prefix>: Prefix
EOF
  exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo '-p <path>: Path where outputs from imputePipe.py are located'
    echo '-s <save-path>: Path where .bgen files are saved'
    shift
    shift
    ;;
    -d|--directory)
    FILEPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -p|--prefix)
    PREFIX="$2"
    shift
    shift
    ;;
    *)    # unknown option
    usage # save it in an array for later
esac
done

# Print
echo "Sorting files in $FILEPATH"

# Get pipeline directory - try multiple methods
if [ -n "${BASH_SOURCE[0]}" ]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    PIPELINE_DIR="$(dirname "$SCRIPT_DIR")"
else
    PIPELINE_DIR="$(pwd)"
fi

# Try to read from settings.json
SETTINGS=""
for test_dir in "$PIPELINE_DIR" "$(pwd)" "/labs/mignot/researchers/smuniz/CASPR2_GWAS/imputePipeline_v2.0.0"; do
    if [ -f "${test_dir}/settings.json" ]; then
        SETTINGS="${test_dir}/settings.json"
        break
    fi
done

if [ -z "$SETTINGS" ]; then
    echo "ERROR: Could not find settings.json"
    exit 1
fi

# Read settings
SLURMDIR=$(jq -r '.folder.SLURM_IMPUTE_LOG' "$SETTINGS")
GWASBYCHR=$(jq -r '.folder.GWAS_BY_CHR' "$SETTINGS")
FILESFOLDER=$(jq -r '.folder.FILESFOLDER' "$SETTINGS" | sed 's|/$||')

# Load module
module load plink2

# Go to path 
cd $FILEPATH

# Find .fam file - try multiple locations
# Note: All chromosomes have the same samples, so we can use the original .fam file
# The cleanup step now keeps .fam files in the original location
FAM_FILE=""
# Try original .fam file first (this is where it stays after cleanup)
if [ -f "$FILESFOLDER/${PREFIX}.fam" ]; then
    FAM_FILE="$FILESFOLDER/${PREFIX}.fam"
# Try in current directory (fallback)
elif [ -f "${PREFIX}.fam" ]; then
    FAM_FILE="${PREFIX}.fam"
# Try chromosome-specific .fam file (legacy support, in case files weren't cleaned up)
elif [ -f "$GWASBYCHR${PREFIX}_CHR${SLURM_ARRAY_TASK_ID}.fam" ]; then
    FAM_FILE="$GWASBYCHR${PREFIX}_CHR${SLURM_ARRAY_TASK_ID}.fam"
else
    echo "ERROR: Could not find .fam file for $PREFIX"
    echo "Searched in:"
    echo "  - $FILESFOLDER/${PREFIX}.fam"
    echo "  - ${PREFIX}.fam"
    echo "  - $GWASBYCHR${PREFIX}_CHR${SLURM_ARRAY_TASK_ID}.fam"
    exit 1
fi

echo "Using .fam file: $FAM_FILE"

# Verify input files exist
INPUT_GEN="CHR${SLURM_ARRAY_TASK_ID}_${PREFIX}.impute.gz"
if [ ! -f "$INPUT_GEN" ]; then
    echo "ERROR: Input file $INPUT_GEN not found!"
    exit 1
fi

# Check input file size
INPUT_SIZE=$(stat -c%s "$INPUT_GEN" 2>/dev/null || echo 0)
if [ "$INPUT_SIZE" -lt 1000000 ]; then
    echo "ERROR: Input file $INPUT_GEN is too small ($INPUT_SIZE bytes). Expected at least 1MB."
    exit 1
fi

echo "Input file: $INPUT_GEN (size: $(numfmt --to=iec-i --suffix=B $INPUT_SIZE))"

# Create sample file (use > to overwrite, not >> to append)
SAMPLE_FILE="CHR${SLURM_ARRAY_TASK_ID}_${PREFIX}.sample"
echo "ID_1 ID_2" > "$SAMPLE_FILE"
echo "0 0" >> "$SAMPLE_FILE"
awk '{
    $1=$2" "$2;
    print $1;
  }' "$FAM_FILE" >> "$SAMPLE_FILE"

# Verify sample file was created correctly
SAMPLE_LINES=$(wc -l < "$SAMPLE_FILE")
if [ "$SAMPLE_LINES" -lt 3 ]; then
    echo "ERROR: Sample file $SAMPLE_FILE has too few lines ($SAMPLE_LINES). Expected at least 3."
    exit 1
fi
echo "Sample file created: $SAMPLE_FILE ($SAMPLE_LINES lines)"

# Convert to bgen using plink2
echo "Starting plink2 conversion for CHR${SLURM_ARRAY_TASK_ID}..."
OUTPUT_PREFIX="CHR${SLURM_ARRAY_TASK_ID}_${PREFIX}"
OUTPUT_BGEN="${OUTPUT_PREFIX}.bgen"

# Redirect plink2 output to a log file to avoid corruption from multiple array jobs
PLINK2_LOG="plink2_CHR${SLURM_ARRAY_TASK_ID}.log"

plink2 \
    --gen "$INPUT_GEN" ref-first \
    --sample "$SAMPLE_FILE" \
    --oxford-single-chr ${SLURM_ARRAY_TASK_ID} \
    --export bgen-1.2 \
    --out "$OUTPUT_PREFIX" > "$PLINK2_LOG" 2>&1

PLINK2_EXIT=$?

# Check the log for any errors
if grep -qi "error\|fatal\|abort\|terminate\|exception" "$PLINK2_LOG"; then
    echo "WARNING: plink2 log contains errors. Last 20 lines:"
    tail -20 "$PLINK2_LOG"
fi

if [ $PLINK2_EXIT -ne 0 ]; then
    echo "ERROR: plink2 failed with exit code $PLINK2_EXIT for CHR${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

# Verify output file was created and has reasonable size
if [ ! -f "$OUTPUT_BGEN" ]; then
    echo "ERROR: Output BGEN file $OUTPUT_BGEN was not created!"
    exit 1
fi

OUTPUT_SIZE=$(stat -c%s "$OUTPUT_BGEN" 2>/dev/null || echo 0)
if [ "$OUTPUT_SIZE" -lt 1000000 ]; then
    echo "ERROR: Output BGEN file $OUTPUT_BGEN is too small ($OUTPUT_SIZE bytes). Expected at least 1MB."
    echo "This suggests the conversion may have failed or only partially completed."
    exit 1
fi

echo "Conversion successful: $OUTPUT_BGEN (size: $(numfmt --to=iec-i --suffix=B $OUTPUT_SIZE))"

