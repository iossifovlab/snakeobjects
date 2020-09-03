#!/bin/bash
#
# the first argument is the name of the project (no '/', '.'  in it)
# the next argument may be '--profile profile_name', where profile_name
# is a folder in $HOME/.config/snakemake/
# profile_name can be slurm_profile and sge_profile
# following the profile you can have the name of a rule or a target
#

name=$1
shift
sbatch --job-name=$name                       \
       --mem=4000                             \
       --error='log/%x-%A_%a.err'             \
       --output='log/%x-%A_%a.out'            \
       run_snake.sh -j1 $*




