#!/bin/bash

if [ -z "$PROJECT_DIR" ]; then
    export PROJECT_DIR=`pwd`
fi

if [ -z "$PIPELINE_DIR" ]; then
    export PIPELINE_DIR=$PROJECT_DIR
fi

echo "WORKING ON PROJECT" $PROJECT_DIR 
echo "WITH PIPELINE" $PIPELINE_DIR

cd $PROJECT_DIR/objLinks

if [ -f "${PROJECT_DIR}/parameters.yaml" ]; then
    default_options=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2`
fi

input=($default_options)
echo ${input[@]}

N=`echo ${input[@]}|wc -w`
profile=""

for n in `seq 0  $N`
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
	echo "profile $profile"
fi

cmd=""

if [ ! -f $profile/config.yaml ]; then
	echo "no config.yaml in the profile $profile"
    exit 1
else
	cmd=`grep -P "^cluster\:" $profile/config.yaml|cut -d':' -f2`
fi

if [ -z "$cmd" ]; then
    echo "no cluster specified in profile config.yaml file"
	exit 1
else
    iippl jobscript.sh >$PROJECT_DIR/jobscript.sh
    echo "$default_options $* && exit 0 || exit 1" >>$PROJECT_DIR/jobscript.sh
    $cmd $PROJECT_DIR/jobscript.sh &
fi

