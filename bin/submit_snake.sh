#!/bin/bash

sbatch --job-name=${PROJECT_DIR}              \
       --mem=4000                             \
       --error='log/%x-%A_%a.err'             \
       --output='log/%x-%A_%a.out' run_snake.sh --profile SLURM.nygc $*




