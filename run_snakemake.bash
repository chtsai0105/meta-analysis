#!/bin/bash

snakemake -p --profile slurm -k --use-conda --use-envmodule --use-singularity --jobs 16 --max-threads 8
