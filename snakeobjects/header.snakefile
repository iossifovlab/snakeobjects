shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, add_targets, get_targets

# CLOUD related reorganization
# add the following ot the main snakefile
#
# if os.envinor["SO_CONTAINER"] == 'yes' and both --default-provider and --default-prefix are in the arguments
#   call download_project_files_from_remote(provider,prefix)

project = Project()

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

