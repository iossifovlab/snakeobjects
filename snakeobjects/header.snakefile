shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, add_targets, get_targets

project = Project()

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

