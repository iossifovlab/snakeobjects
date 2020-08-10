shell.prefix("set -o pipefail; ")
configfile: "OG.json"

rule all:
  input:
    "objLinks/base/o/obj.flag"

rule build_a:
  output:
    "objLinks/base/o/a.txt"
  shell:
    "touch {output}"

rule obj_flag:
  input:
    "objLinks/base/o/a.txt"
  output:
    "objLinks/base/o/obj.flag"
  shell:
    "touch {output}"

