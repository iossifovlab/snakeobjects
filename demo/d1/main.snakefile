shell.prefix("set -o pipefail; ")
configfile: "OG.json"

pipeline = os.environ["PIPELINE_DIR"]

include: pipeline + "/base.snakefile"

rule all_main:
  input:
    "objLinks/base/o/obj.flag"

