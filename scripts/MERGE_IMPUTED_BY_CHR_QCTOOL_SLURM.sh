#!/bin/bash -l
#SBATCH --job-name=MERGE_BY_CHR
#SBATCH --mem-per-cpu=64000
#SBATCH --account=mignot
#SBATCH --time=12:00:00

# Load module
module load qctool/v2.0.1 

# Go to path 
path=/labs/mignot/Oxford_LGI1_Imputed
cd $path 

# Create mergelist 
find ./*.bgen -printf "%f\n" | cut -d '.' -f 1 > mergeList.txt

# Create array 
mapfile -t FILEARRAY < mergeList.txt
LENARRAY=${#FILEARRAY[@]}

# First pass
qctool_v2.0.1 -g ${FILEARRAY[1]}.bgen -s ${FILEARRAY[1]}.sample \
    -g ${FILEARRAY[2]}.bgen -s ${FILEARRAY[2]}.sample \
    -g ${FILEARRAY[3]}.bgen -s ${FILEARRAY[3]}.sample
    -og temp.bgen -os temp.sample 

# # Iterate over for loop
# for IDX in ${FILEARRAY[@]:3:$LENARRAY}; do 

#     # Compute QCTOOOL
#     qctool_v2.0.1 -g temp.bgen -s temp.sample \
#         -merge-in ${FILEARRAY[IDX]}.bgen ${FILEARRAY[IDX]}.sample \
#         -og temp.bgen -os temp.sample
    
# done 
