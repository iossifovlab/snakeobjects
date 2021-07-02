shell.prefix("set -o pipefail; ")
from snakeobjects import Project 
from snakeobjects.snakeUtils import set_project, T, TE, DT, LFS, P, DP, PP, B, add_targets, get_targets
from snakeobjects.remoteProjects import upload_project_files_to_remote
# CLOUD related reorganization
import os,sys
from snakemake import get_argument_parser
parser=get_argument_parser()
args = parser.parse_args(sys.argv)

if ("SO_CONTAINER" in os.environ and
   os.environ["SO_CONTAINER"] == 'yes' and 
   args.default_remote_prefix and 
   args.default_remote_prefix):
	upload_project_files_to_remote(args.default_remote_prefix,
		default_remote_prefix)

project = Project()

set_project(project)

rule so_all_targets:
  input:
    expand("{of}", of=project.get_all_object_flags())

