#!/bin/bash -l
#SBATCH --job-name=CLEAN_FOLDER
#SBATCH --mem-per-cpu=8000
#SBATCH --time=05:00:00
#SBATCH --account=mignot
rm *_diplotype_ordering *_info_by_sample *_summary *_warnings *info
find . -name '*.[[:digit:]][[:digit:]][[:digit:]]' -exec rm -f {} \;
find . -name '*.[[:digit:]][[:digit:]]' -exec rm -f {} \;
find . -name '*.[[:digit:]]' -exec rm -f {} \;
