#!/bin/bash 

#SBATCH --job-name=PLINK_SPLIT_WRAPPER
#SBATCH --output=PLINK_SPLIT_WRAPPER.out
#SBATCH --error=PLINK_SPLIT_WRAPPERsqueue.err
#SBATCH --mem-per-cpu=16000
#SBATCH --account=mignot
#SBATCH --time=1:00:00

command="./bin/plink --bfile "$1" --chr \$SLURM_ARRAY_TASK_ID --make-bed --out "$1"_CHR\$SLURM_ARRAY_TASK_ID"
touch PLINK_SPLIT.sh
chmod 755 PLINK_SPLIT.sh
cat > PLINK_SPLIT.sh <<- EOF
#!/bin/bash -l
#SBATCH --job-name=PLINK_SPLIT_CHR
#SBATCH --mem-per-cpu=4000
#SBATCH --time=01:00:00
#SBATCH --array=1-22
#SBATCH --account=mignot
$command
EOF
sbatch --export=ALL PLINK_SPLIT.sh
