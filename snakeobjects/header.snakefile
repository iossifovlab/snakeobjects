shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, EF, add_targets, get_targets
from snakeobjects.remoteProjects import download_project_files_from_remote

# CLOUD related reorganization
import os,sys
from snakeobjects.cli import get_arg_value
#''' THIS IS BROKEN!!
if os.environ["SO_CONTAINER"] == 'yes'):
  provider = get_arg_value(sargs,'--default-remote-provider')
  bucket = get_arg_value(sargs,'--default-remote-prefix')
  if provider and bucket:
    download_project_files_from_remote(provider,bucket)

if SO_KUBERNETES in os.environ:
  os.environ['SO_PIPELINE'] = "." #or "/source" # ???
  os.environ['SO_PROJECT'] = "./" + bucket 
  fixPP = os.environ['SO_PIPELINE'] + "/fix_permissions.sh"
  if is_file(fixPP):
    os.system("sh fixPP")

#'''	

project = Project()
so_pipeline=project.get_pipeline_directory()

os.environ['PATH']=so_pipeline + ":" + os.environ['PATH']

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

