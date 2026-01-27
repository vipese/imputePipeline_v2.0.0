#!/bin/bash 

command1="rm *_diplotype_ordering *_info_by_sample *_summary *_warnings *info"
command2="find . -name '*.[[:digit:]][[:digit:]][[:digit:]]' -exec rm -f {} \;"
command3="find . -name '*.[[:digit:]][[:digit:]]' -exec rm -f {} \;"
command4="find . -name '*.[[:digit:]]' -exec rm -f {} \;"
touch CLEAN_FOLDER.sh
chmod 755 CLEAN_FOLDER.sh
cat > CLEAN_FOLDER.sh <<- EOF
#!/bin/bash -l
#SBATCH --job-name=CLEAN_IMPUTE
#SBATCH --mem-per-cpu=8000
#SBATCH --time=05:00:00
#SBATCH --account=mignot
$command1
$command2
$command3
$command4
EOF
sbatch --export=ALL CLEAN_FOLDER.sh
