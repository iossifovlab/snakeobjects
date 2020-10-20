shell.prefix("set -o pipefail; ")
import os,re,yaml
from iippl.snakeUtils import set_object_graph, all_obj_dirs, all_obj_types, T, TE, DT, EFS, P, DP, PP, B
from iippl.ObjectGraph import load_object_graph

pipeline = os.environ["PIPELINE_DIR"]

projDir = os.environ["PROJECT_DIR"]
if os.path.exists(projDir + "/parameters.yaml"):
  CF = open(projDir + "/parameters.yaml", 'r')
  if CF:
    config = yaml.safe_load(CF)
  CF.close()
else:
  config = {}

def load_project_params(conf):
  ptn = re.compile('(^\[E\:)(.*)(\])')
  for k,v in conf.items():
    if type(v) != str:
      continue
    m = ptn.match(v)
    if m:
      s = m.span()
      name = m.groups(0)[1]
      n = os.environ[name]
      if n:
        conf[k] = v.replace(v[s[0]:s[1]],n)
      else:
        print('Varianble %s is not defined' % name, file=sys.stderr)
        exit(1) 
  return conf  

if config:
  config = load_project_params(config)

objectGraphFile = ".pipes/OG.json"
if 'objectGraphFile' in config: 
    objectGraphFile = config['objectGraphFile']
OG = load_object_graph(objectGraphFile) 
set_object_graph(OG)

def GP(parameter):
    return config[parameter]


