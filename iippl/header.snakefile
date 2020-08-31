shell.prefix("set -o pipefail; ")
import os
from iippl.snakeUtils import set_object_graph, all_obj_dirs, all_obj_types, T, TE, DT, EFS, P, DP, PP
from iippl.ObjectGraph import load_object_graph

pipeline = os.environ["PIPELINE_DIR"]

if os.path.exists("../parameters.yaml"):
    configfile: "../parameters.yaml"
else:
    config = {}

objectGraphFile = "../OG.json"
if 'objectGraphFile' in config: 
    objectGraphFile = config['objectGraphFile']
OG = load_object_graph(objectGraphFile) 
set_object_graph(OG)

def GP(parameter):
    return config[parameter]


