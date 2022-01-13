#!/bin/bash
# properties = {"type": "single", "rule": "run_snake", "input": {}, "output":{},"params":{}, "log":{}, "wildcards": {}, "local": false, "threads": 1, "resources": {"mem": 64000}, "cluster": {}}
 cd $SO_PROJECT && \
python \
-m snakemake --snakefile $SO_PIPELINE/Snakefile \
--force -j --keep-target-files --keep-remote \
--attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
--nocolor --notemp --no-hooks --use-conda \
--restart-times 3 --max-jobs-per-second 1 --max-status-checks-per-second 1 --local-cores 1 --jobs 800 --latency-wait 90 