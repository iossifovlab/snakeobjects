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

default_args=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2` 

input=($default_args)

N=`echo ${input[@]}|wc -w`
profile=""
for n in `seq 1  $N`
do
	if [[ ${input[$n]} ==  "--profile" ]]; then
		let k=$n+1
		profile=${input[$k]}
	fi
done

if [ -z $profile ]; then
	echo "no profile specified"
	exit 1
else
	echo $profile
fi

cmd=""

if [ ! -f $HOME/.config/snakemake/$profile/config.yaml ]; then
	echo "no config.yaml in the profile $HOME/.config/snakemake/$profile"
	if [ ! -f $profile/config.yaml ]; then
		echo "no config.yaml in the profile $profile"
		exit 1
	else
		cmd=`grep -P "^cluster\:" $profile/config.yaml|cut -d':' -f2`
    fi
else
cmd=`grep -P "^cluster\:" $HOME/.config/snakemake/$profile/config.yaml|cut -d':' -f -2`
fi

if [ -z "$cmd" ]; then
    echo "no cluster specified in profile config.yaml file"
	exit 1
else
	echo "$cmd $PROJECT_DIR/jscript.sh $default_args $* 2>tmp &"
	$cmd jscript.sh $default_args $*  &
fi

