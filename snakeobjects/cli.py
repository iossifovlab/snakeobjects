from argparse import ArgumentParser
from snakeobjects import __version__, Project, ObjectGraph
import importlib.resources as importlib_resources
import os,sys
from importlib.util import spec_from_file_location, module_from_spec


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

    if command in ["prepare","prepareTest"]:
        bldObjGraphPy = proj.get_pipeline_directory() + "/build_object_graph.py"
        if os.path.isfile(bldObjGraphPy):
            spec = spec_from_file_location("build_object_graph", bldObjGraphPy)
            foo = module_from_spec(spec)
            spec.loader.exec_module(foo)
            
            newObjectGraph = ObjectGraph()
            foo.run(proj,newObjectGraph,*args[1:])
        else:
            print(f'ERROR: There is no {bldObjGraphPy}')
            exit(1)
        if command == "prepareTest":
            print("Current graph stats")
            print("+++++++++++++++++++")
            proj.objectGraph.print_stats()
            print("\n")
            print("New graph stats")
            print("+++++++++++++++")
            newObjectGraph.print_stats()
        else:
            proj.objectGraph = newObjectGraph 
            proj.save_object_graph()
            proj.prepare_objects()
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
        print("UPDATING ENVIRONMENT:")
        print("export SO_PROJECT=",proj.directory,sep="") 
        print("export SO_PIPELINE=",proj.get_pipeline_directory(),sep="") 
        print("export PATH=$SO_PIPELINE:$PATH",sep="")
        print("RUNNING:", " ".join(sargs))
        os.environ['SO_PROJECT']  = proj.directory
        os.environ['SO_PIPELINE'] = proj.get_pipeline_directory() 
        os.environ['PATH'] = proj.get_pipeline_directory() + ":" + os.environ['PATH']
        os.execvp('snakemake',sargs)
    elif command == "describe":
        print("Project parameters:")
        for k,v in proj.parameters.items():
            print(f"\t{k}: {v}")
        proj.objectGraph.print_stats()
    elif command == "version":
        return __version__
    elif command in ["help", "-h", "--help"]:
        print("Available commands are: \033[1m prepare, prepareTest, prepareObjects, run, describe, and version.\033[0m\n")
        print("Typical sequence of commands is descripe, prepareTest, prepare, run:\n")
        print("\033[1mdescrbe\033[0m Prints a basic information about the project and the pipline that are used;\n")
        print("\033[1mprepareTest\033[0m Uses the build_object_graph.py in the pippeline directory to create a new object graph and prints statistics of the current and the new object graph. The project is not modified at all;\n")
        print("\033[1mprepare\033[0m First, uses the build_object_graph.py in the pippeline directory to create an object graph for the snakeobjects project. Then prepares the main snakefile, and the directories and symbolic links for all object in the object graph;\n")
        print("\033[1mprepareObjects\033[0m Prepares the main snakefile, and the directories and symbolic links for all object in the object graph;\n")
        print("\033[1mrun [<arguments to snakemake>]\033[0m Creates targets for object in the object graph by running snakemake. The <arguments to snakemake>  determine which targets will be created and what resources will be used;\n")
        print("\033[1mversion\033[0m This prints the version;\n")
        print("\033[1mhelp\033[0m This shows help;\n")
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
