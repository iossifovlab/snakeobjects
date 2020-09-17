shell.prefix("set -o pipefail; ")
import os
from iippl.snakeUtils import set_object_graph, all_obj_dirs, all_obj_types, T, TE, DT, EFS, P, DP, PP, B
from iippl.ObjectGraph import load_object_graph

pipeline = os.environ["PIPELINE_DIR"]

if os.path.exists("../parameters.yaml"):
    configfile: "../parameters.yaml"
else:
    config = {}

objectGraphFile = ".pipes/OG.json"
if 'objectGraphFile' in config: 
    objectGraphFile = config['objectGraphFile']
OG = load_object_graph(objectGraphFile) 
set_object_graph(OG)

def GP(parameter):
    return config[parameter]


rule all_main:
  input:
    expand("{od}/obj.flag", od=all_obj_dirs())

include: pipeline + "/base.snakefile"

rule all_base:
  input:
    expand("base/{o}/obj.flag", o=OG.names("base"))

include: pipeline + "/level1.snakefile"

rule all_level1:
  input:
    expand("level1/{o}/obj.flag", o=OG.names("level1"))

include: pipeline + "/level2.snakefile"

rule all_level2:
  input:
    expand("level2/{o}/obj.flag", o=OG.names("level2"))

include: pipeline + "/level3.snakefile"

rule all_level3:
  input:
    expand("level3/{o}/obj.flag", o=OG.names("level3"))

