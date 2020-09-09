#!/bin/bash
# properties = {"rule":"run_snake", "threads":1, "mem_free":4000M, "cluster": {}}
 /gpfs/commons/home/byamrom/anaconda3/envs/pipesDev/bin/python3.8 -m \
 snakemake --snakefile $PROJECT_DIR/objLinks/main.snakefile  $*      \
 && exit 0 || exit 1
