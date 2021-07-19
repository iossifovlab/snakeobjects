shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, EF, add_targets, get_targets
from snakeobjects.remoteProjects import download_project_files_from_remote

# CLOUD related reorganization
import os,sys

if ("SO_CONTAINER" in os.environ and 
  "SO_KUBERNETES" in os.environ and
  os.environ["SO_CONTAINER"] == 'yes'):
    provider, bucket = os.environ["SO_KUBERNETES].split(":")
    download_project_files_from_remote(provider,bucket)
    os.environ['SO_PIPELINE'] = "."
    os.environ['SO_PROJECT'] = "./" + bucket 
    fixPP = os.environ['SO_PIPELINE'] + "/fix_permissions.sh"
    if os.path.exists(fixPP):
      os.system("sh fixPP")

project = Project()
so_pipeline=project.get_pipeline_directory()

os.environ['PATH']=so_pipeline + ":" + os.environ['PATH']

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

