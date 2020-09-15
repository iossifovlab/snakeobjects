#!/bin/bash
# properties = {"type": "single", "rule": "run_snake", "input": {}, "output":{},"params":{}, "log":{}, "wildcards": {}, "local": false, "threads": 1, "resources": {"mem": 8000}, "cluster": {}}
 cd $PROJECT_DIR/objLinks && \
python \
-m snakemake --snakefile $PROJECT_DIR/objLinks/.pipes/main.snakefile \
--force -j --keep-target-files --keep-remote \
--latency-wait 40 \
--attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
--nocolor --notemp --no-hooks --nolock \
--mode 2  \
