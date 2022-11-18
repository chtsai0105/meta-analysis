#!/bin/bash

snakemake -p --profile slurm --use-conda --use-envmodule --jobs 6 --max-threads 20