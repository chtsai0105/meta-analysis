#!/bin/bash

for round in {1..5}; do
    snakemake -p --profile slurm -k --use-conda --use-envmodule --use-singularity --jobs 16 --max-threads 8 --batch all=${round}/5
done
