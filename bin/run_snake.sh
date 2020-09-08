#!/bin/bash

if [ -z "$PROJECT_DIR" ]; then
    export PROJECT_DIR=`pwd`
fi

if [ -z "$PIPELINE_DIR" ]; then
    export PIPELINE_DIR=$PROJECT_DIR
fi

echo "WORKING ON PROJECT" $PROJECT_DIR 
echo "WITH PIPELINE" $PIPELINE_DIR

mkdir -p $PROJECT_DIR/objLinks
cd $PROJECT_DIR/objLinks

if [ ! -f header.snakefile ]; then

    iippl header.snakefile > $PROJECT_DIR/header.snakefile
fi

default_opetions=`grep default_snakemake_args parameters.yaml |cut -d'=' -f2`
if [ -z "$default_options" ]; then
snakemake --snakefile main.snakefile  $*
else
snakemake --snakefile main.snakefile  $default_options $*
fi

# --profile $BIN_DIR/SLURM.nygc $*
