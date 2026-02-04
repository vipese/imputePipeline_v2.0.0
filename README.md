# Imputation_pipeline_v2

## Introduction
This pipeline attempts to update and improve the [1000 Genome Imputation Pipeline](https://github.com/Mignot-Lab/imputePipeline) designed by [@adiamb](https://github.com/adiamb). This version automates the process and provides a final output of imputed PED binary (_.bed_) files using a settings-logic. 

This pipeline is SLURM-dependant, and therefore may only be used with Stanford SCG cloud-computing services. 

## Description
The pipeline performs the following steps:
1. Takes base PLINK format _.bed_ files, splits them into chromosomes (1 to 22), phases them using [SHAPEIT](https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html), and subsequently imputes each chromosome individually to [1000 Genomes Phase I](https://www.internationalgenome.org/) in segments of 1 mb. (Unmodified  from [Original Code](https://github.com/Mignot-Lab/imputePipeline)).
2. Cleans up SLURM and SHAPEIT _.log_ files (Added to Original Code).
3. Concatenates the imputed segments from step 1) (Modified from Original Code).
4. Sorts and converts to _.bgen_ format (Modified from Original Code).

## Requirements
### Packages
The pipeline requires the following modules: <br>
``plink``  <br>
``qctool/v2.0.1``

### Settings
To operate the pipeline, all the arguments must be introduced through the `settings.json` file. Particularly, you must keep in mind:
1. To write the __total path__, and _not_ the relative path.
2. To write a slash (/) at the end of each path. 

#### Folder Settings

| Setting | Description |
|---------|-------------|
| `FILESFOLDER` | Folder containing the pipeline (repository) |
| `SOURCE_DATA` | Folder containing non-imputed PLINK files (used by `submit_batch_imputation.sh`) |
| `GWAS_BY_CHR` | Folder for intermediate chromosome-split files |
| `SLURM_IMPUTE_LOG` | Folder for SLURM output and error logs |
| `SHAPEIT_IMPUTE_LOG` | Folder for SHAPEIT logs |
| `BIN_FOLDER` | Folder for imputed output files |

For more information regarding each of the elements in the `settings.json` file, please refer to the _info_ element within the file.

#### Pipeline Options (in main.sh)

| Variable | Default | Description |
|----------|---------|-------------|
| `CLEANUP` | `true` | Remove intermediate files after pipeline completion |
| `CLEANUP_CHR_FILES` | `false` | Remove per-chromosome files from output folder after merge |

## Usage

### Single File Imputation
To run the pipeline for a single PLINK file:
1. Clone the repository.
2. Copy the _.bed_, _.bim_, and _.fam_ files inside the directory.
3. Fill up the `settings.json` file (set `prefix` to your file prefix).
4. Run from the repository's directory: `sbatch main.sh`

### Batch Imputation
To impute multiple PLINK files automatically, use the batch submission script:

```bash
# Preview what would be submitted (no actual jobs)
./submit_batch_imputation.sh --dry-run

# Submit all unimputed files
./submit_batch_imputation.sh

# Test with first 3 files only
./submit_batch_imputation.sh --test

# Combine flags
./submit_batch_imputation.sh --dry-run --test
```

The batch script:
1. Reads `SOURCE_DATA` from `settings.json` to find non-imputed PLINK files
2. Checks which files have already been imputed (exist in `BIN_FOLDER`)
3. Submits a separate SLURM job for each unimputed file
4. Each job runs in an isolated working directory to prevent conflicts
5. Results are copied to `BIN_FOLDER` upon completion

Logs are saved to `BATCH_IMPUTE_LOGS/` with timestamps to prevent overwrites.

## Issues 
There are certain issues that must be taken into consideration prior to utilizing the pipeline:
1. Low memory: the pipeline can be computationally heavy, and may not perform correctly with large datasets. Allow always ~50 GB of free space.
2. The pipeline will fail if the files to be imputed contain duplicated variants (by position) or IDs.
3. As of yet, the pipeline does not convert to _.bed_ files. This is because it converts to _.bgen_ and to allows merging it with other databases.
4. Pipeline may not function if the file to impute has multi-allelic variants.
5. There is an issue when merging by chromosomes using QCTOOLS. When it converts to _.fam_ it messes up the IIDs.

## Future changes
The following changes are recommended to improve the pipeline:
1. Add headers to files.
2. Improve `main.sh`
3. Cleaning may be put into a single bash file in order to tidy up `main.sh`
4. Include QC (remove duplicated variants and IIDs).
5. Finalize *utils*, to merge by CHR and between datasets. Maybe include it in main and add `settings.json` option.
6. Add more verbose.

