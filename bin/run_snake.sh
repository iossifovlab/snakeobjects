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

if [ ! -f main.snakefile ]; then
    iippl main.snakefile > main.snakefile
fi

snakemake --snakefile main.snakefile  $*

# --profile $BIN_DIR/SLURM.nygc $*
