shell.prefix("set -o pipefail; ")
configfile: "OG.json"

pipeline = os.environ["PIPELINE_DIR"]

include: pipeline + "/base.snakefile"


RNASEQs = sorted([x[:-7] for x in config if x.endswith(".rnaSeq")])
SEQLIBs = sorted([x[:-7] for x in config if x.endswith(".seqLib")])

rule all_main:
  input:
    "objLinks/base/o/obj.flag"

