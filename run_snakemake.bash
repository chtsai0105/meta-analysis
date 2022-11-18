#!/bin/bash

snakemake -p --profile slurm --use-envmodule --jobs 6 --max-threads 20