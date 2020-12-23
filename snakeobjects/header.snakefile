shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B

project = Project()

set_project(project)

rule all_main:
  input:
    expand("{of}", of=project.get_all_object_flags())

