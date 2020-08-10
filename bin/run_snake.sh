#!/bin/bash

mkdir -p $PROJECT_DIR/objLinks
cd $PROJECT_DIR/objLinks
snakemake --snakefile <(iippl main.snakefile)  $*

# --profile $BIN_DIR/SLURM.nygc $*
