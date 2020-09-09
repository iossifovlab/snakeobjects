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

if [ ! -f header.snakefile ]; then

    iippl header.snakefile > $PROJECT_DIR/header.snakefile
fi

cmd=`grep -P "^cluster" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2`
echo $cmd
default_args=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2` 

if [ -z "$cmd" ]; then
    echo "run_snake.sh -j1 $*"
	run_snake.sh -j1 $* 2> tmp &
else
	echo "$cmd $PROJECT_DIR/jscript.sh $default_args $* 2>tmp &"
	$cmd $PROJECT_DIR/jscript.sh $default_args $* 2> tmp &
fi




