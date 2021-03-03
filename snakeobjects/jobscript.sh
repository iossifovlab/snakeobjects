#!/bin/bash
# properties = {"type": "single", "rule": "run_snake", "input": {}, "output":{},"params":{}, "log":{}, "wildcards": {}, "local": false, "threads": 1, "resources": {"mem": 8000}, "cluster": {}}
 cd $SO_PROJECT/objects && \
python \
-m snakemake --snakefile $SO_PROJECT/objects/.snakeobjects/main.snakefile \
--force -j --keep-target-files --keep-remote \
--attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
--nocolor --notemp --no-hooks --use-conda \
