#!/bin/bash

snakemake -p --profile slurm --use-conda --use-envmodule --jobs 16 --max-threads 8