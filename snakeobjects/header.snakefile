shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, add_targets, get_targets

# CLOUD related reorganization
# add the following ot the main snakefile
#
# if os.envinor["SO_USE_REMOTE"] == 'yes':
#   get the value of  --default-provider and --default-prefix from snakemake somehow.
#   call download_project_files_from_remote(provider,prefix)


project = Project()

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

