shell.prefix("set -o pipefail; ")
import os,re,yaml
from snakeobjects.snakeUtils import set_object_graph, all_obj_dirs, all_obj_types, T, TE, DT, EFS, P, DP, PP, B
from snakeobjects.ObjectGraph import load_object_graph, load_project_params

pipeline = os.environ["PIPELINE_DIR"]

projDir = os.environ["PROJECT_DIR"]
if os.path.exists(projDir + "/parameters.yaml"):
  config = load_project_params(projDir + "/parameters.yaml")
else:
  config = {}

objectGraphFile = ".snakeobjects/OG.json"
if 'objectGraphFile' in config: 
    objectGraphFile = config['objectGraphFile']
OG = load_object_graph(objectGraphFile) 
set_object_graph(OG)

def GP(parameter):
    return config[parameter]


