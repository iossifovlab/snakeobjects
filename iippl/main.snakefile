shell.prefix("set -o pipefail; ")
configfile: "../OG.json"
pipeline = os.environ["PIPELINE_DIR"]
import os

from iippl.snakeUtils import set_config, all_obj_ids, all_obj_dirs, all_obj_types, T, DT, EFS, P, DP, PP

set_config(config)

GP = config['parameters']

for t in all_obj_types():
    tsfn = pipeline + "/" + t + ".snakefile"
    if not os.path.isfile(tsfn):
        print("Creating a starter", tsfn)
        F = open(tsfn,"w")
        print("rule %s_obj:\n"
              "  input:\n"
              "    DT('obj.flag')\n"
              "  output: \n"
              "    touch(T('obj.flag'))" % (t),file=F)
        F.close()
    include: tsfn

rule all_main:
  input:
    expand("{od}/obj.flag", od=all_obj_dirs())

