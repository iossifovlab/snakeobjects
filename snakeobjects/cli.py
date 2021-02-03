from argparse import ArgumentParser
from snakeobjects import __version__, Project, ObjectGraph, load_object_graph, graph
import importlib.resources as importlib_resources
import os,sys,yaml
from importlib.util import spec_from_file_location, module_from_spec

def load_yaml(file_name):

    CF = open(file_name, 'r')
    config = yaml.safe_load(CF)
    CF.close()

    return config  

helpData = {
    "version": "prints the version",

    "help":    '''sobjects help [command]

Shows description of the command line interface. 
Without an argument, prints help for all the commands. If one argument 
is given, it should be one of the available commands and help for the given 
command is shown.''',

    "describe": '''Prints a basic information about the project and the pipeline
that are used''',

    "prepareTest": '''sobjects prepareTest [<arguments to build_object_graph.py>]

Uses the build_object_graph.py in the pipeline directory to create a new object 
graph and prints statistics of the current and the new object graph. The project 
is not modified at all.''',

    "prepare": '''sobjects prepare [<arguments to build_object_graph.py>]

First, uses the build_object_graph.py in the pipeline directory 
to create an object graph for the snakeobjects project. Then prepares the main 
snakefile, and the directories and symbolic links for all object in the object 
graph''',

    "prepareObjects": '''sobjects prepareObjects

Prepares the main snakefile, and the directories and symbolic links for all 
object in the object graph''',

    "run": '''sobjects run [<arguments to snakemake>]

Creates targets for object in the object graph by running snakemake. The 
<arguments to snakemake> determine which targets will be created and what resources 
will be used''',

    "submit": '''sobjects submit [<arguments to snakemake>]

Creates targets for object in the object graph by running snakemake with profile specified in default_snakemake_args directive of so_project.yaml. The 
<arguments to snakemake> determine which targets will be created and what resources 
will be used''',

    "graph": '''sobjects graph [<width> <penwidth> <arrowsize>] [<objectGraph json file>]

Prints to stdout dot file for project ObjectGraph. Defaults values for width, penwidth, and arrowsize are 0.05, 0.1, 0.1. If fourth argument is provided, the dot file for the object graph in json file will be created'''

}


def cli(args=None):
    if not args:
        args = sys.argv[1:]
    command = args[0] 


    if command == "jobscript.sh":
        print(importlib_resources.read_text(__package__,'jobscript.sh'),end='')
        return

    if command == "version":
        print(__version__)
        return
    elif command in ["help", "-h", "--help"]:
        print("snakeobject %s\n" % (__version__))
        if len(args) == 1:
            print("Available commands are:\n\t", "\n\t".join(helpData.keys()),"\n",sep="")
        # print("Typical sequence of commands is descripe, prepareTest, prepare, run:\n")
            for cmd,hs in helpData.items():
                print(cmd)
                print('#' * len(cmd))
                print(hs)
                print()
        elif len(args) == 2:
            hCmd = args[1]
            if hCmd in helpData:
                print(helpData[hCmd])
            else:
                print("The command", hCmd, "is unknown")
                return 1
        else:
            print("Help accepts at most one argument.")
            return 1
        return

    proj = Project()
    print("# WORKING ON PROJECT", proj.directory)
    print("# WITH PIPELINE", proj.get_pipeline_directory())

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
    elif command == "submit":
        sargs = []
        if "default_snakemake_args" in proj.parameters:
            sargs += proj.parameters["default_snakemake_args"].split()
            sargs += args[1:]
            print("sargs:", sargs)
        else:
            raise ProjectException("No profile specified")
        profile=sargs[sargs.index('--profile')+1]
        if not os.path.exists(profile): 
            raise ProjectException("Profile not found %s" % profile)
        if not os.path.exists(profile+"/config.yaml"): 
            raise ProjectException("No config.yaml in %s" % profile)
        pr_config = load_yaml(profile+"/config.yaml")
        if not "cluster" in pr_config: 
            ProjectException("cluster in not specified in %s" % profile+"/config.yaml")
        cmd=pr_config["cluster"]
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
        if os.system('sobjects jobscript.sh >$SO_PROJECT/objects/.snakeobjects/jobscript.sh'):
            raise ProjectException("sobjects jobscript.sh failed")
        os.system("%s/%s" % (profile,cmd)+ " $SO_PROJECT/objects/.snakeobjects/jobscript.sh")
        #os.execvp('python', [profile + "/" +cmd, "$SO_PROJECT/objects/.snakeobjects/jobscript.sh"])        
    elif command == "describe":
        print("Project parameters:")
        for k,v in proj.parameters.items():
            print(f"\t{k}: {v}")
        proj.objectGraph.print_stats()
    elif command == "graph":
        print(args, file=sys.stderr)
        if len(args) == 1:
            graph.plotGraph(proj.objectGraph)
        elif len(args) == 4:
            w, pw, ar = [float(a) for a in args[1:]]
            graph.plotGraph(proj.objectGraph, width=w, penwidth=pw, arrowsize=ar)
        elif len(args) == 5:
            w, pw, ar = [float(a) for a in args[1:-1]]
            graph.plotGraph(load_object_graph(args[-1]), width=w, penwidth=pw, arrowsize=ar)
    else:
        print("Don't know the command:", command)
        return 1

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
