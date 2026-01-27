#!/bin/bash -l
#SBATCH --job-name=QCTOOL_CONVERT_BED
#SBATCH --output=QCTOOL_CONVERT_BED.out
#SBATCH --error=QCTOOL_CONVERT_BED.err
#SBATCH --mem-per-cpu=64000
#SBATCH --time=120:00:00
#SBATCH --account=mignot
#SBATCH --array=1-22

# Parse arguments
PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-p <path>]

-p <path>: Path where .bgen files are located
-s <save-path>: Path where PED binary files will be saved
EOF
  exit 1
}

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    echo '-p <path>: Path where .bgen files are located'
    echo '-s <save-path>: Path where PED binary files will be saved'
    shift
    shift
    ;;
    -p|--path)
    PATH="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--save-path)
    SAVEPATH="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    usage() # save it in an array for later
    ;;
esac

# Read settings
SLURMDIR=$(jq -r '.folder.SLURM_IMPUTE_LOG' settings.json)

# Load Module
module load qctool/v2.0.1

# Go to path
cd $PATH

# Create directory if does not exist 
[ -d $SAVEPATH ] || mkdir $SAVEPATH   

# Convert to binary PED
qctool_v2.0.1 \
-filetype bgen \
-precision 2 \
-threads 16 \
-threshold 0.7 \
-assume-chromosome $SLURM_ARRAY_TASK_ID \
-g CHR"$SLURM_ARRAY_TASK_ID"_LGI1_QC.bgen \
-s LGI1_Oxford.sample \
-og  CHR"$SLURM_ARRAY_TASK_ID"_LGI1_QC \
-ofiletype binary_ped

# FAM file does not keep IID, use sample as FAM
awk 'NR>=122 {print $0}' \
CHR"$SLURM_ARRAY_TASK_ID"_LGI1_QC.sample > \
CHR"$SLURM_ARRAY_TASK_ID"_LGI1_QC.fam

# Move files 
mv CHR"$SLURM_ARRAY_TASK_ID"_LGI1_QC* $SAVEPATH
mv QCTOOL_CONVERT_BED* $SLURMDIR
