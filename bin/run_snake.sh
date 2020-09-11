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
mkdir -p log

default_options=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2`

if [ -z "$default_options" ]; then

echo "snakemake --snakefile .pipes/main.snakefile  $*"
snakemake --snakefile .pipes/main.snakefile  $*
else

echo "snakemake --snakefile .pipes/main.snakefile  $default_options $*"
snakemake --snakefile .pipes/main.snakefile  $default_options $*
fi

# --profile $BIN_DIR/SLURM.nygc $*
