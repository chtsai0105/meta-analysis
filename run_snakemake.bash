#!/bin/bash

## Launch all jobs at once.
snakemake -p --profile slurm -k --use-conda --use-envmodule --use-singularity --jobs 16 --max-threads 8
## Please note that the --keep-going (-k) option is enabled to avoid the entire workflow being stopped by a single job failure.
## Make sure to check the failure jobs afterward.


## For batch run. If you have tons of samples and jobs (more than 10000 samples), it is recommended to use batch run.
## However, it would have some trouble at the final batch if some of the jobs failed. Please use in cautious.
# for round in {1..5}; do
#     snakemake -p --profile slurm -k --use-conda --use-envmodule --use-singularity --jobs 16 --max-threads 8 --batch all=${round}/5
# done
