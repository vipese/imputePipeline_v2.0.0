#!/bin/bash 

#SBATCH --job-name=CAT_IMPUTEFILES
#SBATCH --output=CAT_IMPUTEFILES.out
#SBATCH --error=CAT_IMPUTEFILESsqueue.err
#SBATCH --mem-per-cpu=16000
#SBATCH --array=1-22
#SBATCH --account=mignot
#SBATCH --time=12:00:00


PROGNAME=$0

usage() {
  cat << EOF >&2
Usage: $PROGNAME [-p <path>]

One or more of the following flags is missing:

-p <path>: Path where outputs from imputePipe.py are located
-s <save-path>: Path where .gz concat files will be located
-c <scripts>: Path where the scripts are located
-p <prefix>: Prefix of the file used for imputation
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
    echo '-s <save-path>:  Path where .gz concat files will be located'
    shift
    shift
    ;;
    -d|--directory)
    FILESPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -s|--save-dir)
    SAVEPATH="$2"
    shift # past argument
    shift # past value
    ;;
    -c|--scripts)
    SCRIPTPATH="$2"
    shift
    shift
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

# Read settings and pwd
SLURMDIR=$(jq -r '.folder.SLURM_IMPUTE_LOG' settings.json)
PWD=$(pwd)

# Go to path 
cd $FILESPATH

# Run Concatenation
python ${SCRIPTPATH}CONCAT_IMPUTE.py -F CHR"${SLURM_ARRAY_TASK_ID}"_${PREFIX}

# Move file 
[ -d $SAVEPATH ] || mkdir $SAVEPATH   # Create directory if does not exist 
mv CHR"${SLURM_ARRAY_TASK_ID}"_$PREFIX*.gz $SAVEPATH
mv ${pwd}/CAT_IMPUTEFILES* $SLURMDIR
