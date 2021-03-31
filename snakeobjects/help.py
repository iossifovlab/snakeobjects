import importlib.resources as importlib_resources
import os,sys
from importlib.util import spec_from_file_location, module_from_spec

helpData = {
    "version": "prints the version",

    "help":    '''sobjects help [command]

Shows description of the command line interface. 
Without an argument, prints help for all the commands. If one argument 
is given, it should be one of the available commands and help for the given 
command is shown.''',

    "describe": '''Prints basic information about the project and the pipeline
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

Prepares the main.snakefile, and the directories and symbolic links for all 
objects in the object graph''',

    "run": '''sobjects run [<arguments to snakemake>]

Creates targets for objects in the object graph by running snakemake. The 
<arguments to snakemake> determine which targets will be created and what resources 
will be used''',

    "submit": '''sobjects submit [<arguments to snakemake>]

Creates targets for objects in the object graph by running snakemake with
profile specified in default_snakemake_args directive of so_project.yaml. The
<arguments to snakemake> determine which targets will be created and what
resources will be used''',

    "graph": '''sobject graph [-w width] [-p penwidth] [-a arrowsize] [-l legend] [-o out] [-i id] [-s shape]

optional arguments:
  -h, --help            show this help message and exit
  -w width, --width width
                        width of node, default is 0.75
  -p penwidth, --penwidth penwidth
                        thickness of edges, default is 1.0
  -a arrowsize, --arrowsize arrowsize
                        multiplicative scale factor for arrowheads, default is 1.0
  -l legend, --legend legend
                        Name of the output legend file, default is no legend
  -o out, --out out     name of the output file, default is stdout
  -t text, --text text  place text in nodes: [|oId|oType:oId|params], default no text
  -s shape, --shape shape
                        shape of the node, default is circle, for all shape names see https://www.graphviz.org/doc/info/shapes.html'''

}

def help(args=None):
    if not args:
        args = sys.argv[1:]
    command = args[0]

    if command in ["help", "-h", "--help"]:
        print("Snakeobjects %s\n" % (__version__))
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
    return

if __name__ == '__main__':
    import sys
    help(sys.argv[1:])
