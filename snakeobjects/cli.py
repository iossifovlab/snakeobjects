from argparse import ArgumentParser
from snakeobjects import __version__, Project
import importlib.resources as importlib_resources
import os,sys

def cli(args=None):
    if not args:
        args = sys.argv[1:]
    command = args[0] 


    if command == "jobscript.sh":
        print(importlib_resources.read_text(__package__,'jobscript.sh'),end='')
        return


    proj = Project()
    print("WORKING ON PROJECT", proj.directory)
    print("WITH PIPELINE", proj.get_pipeline_directory())
    os.environ['SO_PROJECT']  = proj.directory
    os.environ['SO_PIPELINE'] = proj.get_pipeline_directory() 

    if command in ["prepare","prepareTest"]:
        bldObjGraphPy = proj.get_pipeline_directory() + "/build_object_graph.py"
        if os.path.isfile(bldObjGraphPy):
            bargs = [bldObjGraphPy] + args[1:] + [command]
            print("RUNNING: ", " ".join(bargs))
            os.execvp(bldObjGraphPy,bargs)
        else:
            print(f'ERROR: There is no {bldObjGraphPy}')
            exit(1)
    elif command in ["prepareObjects"]:
        proj.prepare_objects()
    elif command == "run":
        sargs = ['snakemake',
                        '-s', proj.directory + '/objects/.snakeobjects/main.snakefile', 
                        '-d', proj.directory + '/objects'] 
        if "default_snakemake_args" in proj.parameters:
            sargs += proj.parameters["default_snakemake_args"].split()
        sargs += args[1:]
        # os.chdir(proj.directory + '/objects')
        print("RUNNING:", " ".join(sargs))
        os.execvp('snakemake',sargs)
    elif command == "describe":
        print("Project parameters:")
        for k,v in proj.parameters.items():
            print(f"\t{k}: {v}")
        proj.objectGraph.print_stats()
    else:
        print("Don't know the command:", command)
        return 1

    # No return value means no error.
    # Return a value of 1 or higher to signify an error.
    # See https://docs.python.org/3/library/sys.html#sys.exit

'''
if [ -z "$PROJECT_DIR" ]; then
    export PROJECT_DIR=`pwd`
fi

if [ -z "$PIPELINE_DIR" ]; then
    export PIPELINE_DIR=$PROJECT_DIR
fi

echo "WORKING ON PROJECT" $PROJECT_DIR 
echo "WITH PIPELINE" $PIPELINE_DIR

cd $PROJECT_DIR/objects

if [ -f "${PROJECT_DIR}/parameters.yaml" ]; then
    default_options=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2`
fi
echo "snakemake --snakefile .snakeobjects/main.snakefile  $default_options $*"
snakemake --snakefile .snakeobjects/main.snakefile  $default_options $*
'''

if __name__ == '__main__':
    import sys
    print("BBBBBBB")
    cli(sys.argv[1:])
