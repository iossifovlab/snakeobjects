shell.prefix("set -o pipefail; ")
configfile: "OG.json"
pipeline = os.environ["PIPELINE_DIR"]
from iippl.snakeUtils import set_config, T, all_obj_ids, all_obj_dirs, all_obj_types

set_config(config)

for t in all_obj_types():
    include: pipeline + "/" + t + ".snakefile"

rule all_main:
  input:
    expand("{od}/obj.flag", od=all_obj_dirs())

