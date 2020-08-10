#!/bin/bash


mkdir -p $PROJECT_DIR/objLinks
cd $PROJECT_DIR/objLinks
snakemake --snakefile $PIPELINE_DIR/main.snakefile  $*

# --profile $BIN_DIR/SLURM.nygc $*
