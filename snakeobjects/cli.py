from argparse import ArgumentParser
from snakeobjects import __version__
import importlib.resources as importlib_resources
import os,sys

def cli(args=None):
    if not args:
        args = sys.argv[1:]
    command = args[0] 

    if command == "header.snakefile":
        print(importlib_resources.read_text(__package__,'header.snakefile'),end='')
    elif command == "jobscript.sh":
        print(importlib_resources.read_text(__package__,'jobscript.sh'),end='')
    elif command == "run":
        if 'PROJECT_DIR' not in os.environ:
            os.environ['PROJECT_DIR'] = os.cwd()
        if 'PIPELINE_DIR' not in os.environ:
            os.environ['PIPELINE_DIR'] = os.environ['PROJECT_DIR'] 
        print("WORKING ON PROJECT",os.environ['PROJECT_DIR'])
        print("WITH PIPELINE",os.environ['PIPELINE_DIR'])
        sargs = ['snakemake','--snakefile', '.snakeobjects/main.snakefile'] + args[1:]
        os.chdir(os.environ['PROJECT_DIR'] + '/objects')
        ## if [ -f "${PROJECT_DIR}/parameters.yaml" ]; then default_options=`grep -P "^default_snakemake_args" ${PROJECT_DIR}/parameters.yaml |cut -d':' -f2` fi
        print(" ".join(sargs))
        os.execvp('snakemake',sargs)
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
