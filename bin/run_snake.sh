#!/bin/bash

snakemake --snakefile $PIPELINE_DIR/main.snakefile  $*

# --profile $BIN_DIR/SLURM.nygc $*
