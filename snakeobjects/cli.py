import os,sys
from snakeobjects.remoteProjects import upload_project_files_to_remote

helpData = {
    "version": "prints the version",

    "help":    '''sobjects help [command]

Shows description of the command line interface. 
Without an argument, prints help for all the commands. If one argument 
is given, it should be one of the available commands and help for the given 
command is shown.''',

    "describe": '''Prints basic information about the project and the pipeline
that are used.''',

    "prepareTest": '''sobjects prepareTest [<arguments to build_object_graph.py>]

Uses the build_object_graph.py in the pipeline directory to create a new object 
graph and prints statistics of the current and the new object graph. The project 
is not modified at all.''',

    "buildObjectGraph": '''sobjects buildObjectGraph [<arguments to build_object_graph.py>]

Uses the build_object_graph.py in the pipeline directory to create a new object 
graph and stores it on the file OG.json in the project directory.''',

    "createMain": '''sobjects createMain

Uses the objectGraph OG.json to write in the pipeline directory the Snakefile.''',

    "createSymbolicLinks": '''sobjects createSymbolicLinks

Uses the objectGraph OG.json to create object directories for objects that have symbolic 
links in object's parameters.''',

    "prepare": '''sobjects prepare [<arguments to build_object_graph.py>]

First, uses the build_object_graph.py in the pipeline directory 
to create an object graph for the snakeobjects project. Then prepares the main 
snakefile (Snakefile in the pipeline directory), and the directories and symbolic links 
for all object in the object graph that have symbolic links in their parameters.''',

    "prepareObjects": '''sobjects prepareObjects

Create directories for objects in the object graph with symbolic links in their parameters.''',

    "run": '''sobjects run [<arguments to snakemake>]

Creates targets for objects in the object graph by running snakemake. The 
<arguments to snakemake> determine which targets will be created and what resources 
will be used.''',

    "submit": '''sobjects submit [<arguments to snakemake>]

Creates targets for objects in the object graph by running snakemake with
profile specified in default_snakemake_args directive of so_project.yaml. The
<arguments to snakemake> determine which targets will be created and what
resources will be used.''',

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

def load_yaml(file_name):

    CF = open(file_name, 'r')
    config = yaml.safe_load(CF)
    CF.close()

    return config  

def get_arg_value(args,arg):
            try:
                i = args.index(arg)
            except ValueError:
                return None 
            if i+1 >= len(args):    
                return None
            return args[i+1]
        
def cli(args=None):
    if not args:
        args = sys.argv[1:]
    command = args[0] 


    if command == "jobscript.sh":
        print(importlib_resources.read_text(__package__,'jobscript.sh'),end='')
        return

    if "READTHEDOCS" in os.environ:
        from _version import get_versions
        __version__ = get_versions()['version']
    else:
        from snakeobjects import __version__


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
                return 1
        else:
            print("Help accepts at most one argument.")
            return 1
        return


    if command == "version":
        print(__version__)
        return

    from snakeobjects import Project, ObjectGraph, load_object_graph, graph
    import importlib.resources as importlib_resources
    import yaml
    from importlib.util import spec_from_file_location, module_from_spec

    proj = Project()
    print("# WORKING ON PROJECT", proj.directory)
    print("# WITH PIPELINE", proj.get_pipeline_directory())
    def buildObjectGraph():
        bldObjGraphPy = proj.get_pipeline_directory() + "/build_object_graph.py"
        if os.path.isfile(bldObjGraphPy):
            spec = spec_from_file_location("build_object_graph", bldObjGraphPy)
            foo = module_from_spec(spec)
            spec.loader.exec_module(foo)
            newObjectGraph = ObjectGraph()
            foo.run(proj,newObjectGraph,*args[1:])
            return newObjectGraph 
        else:
            print(f'ERROR: There is no {bldObjGraphPy}')
            exit(1)
    if command == "buildObjectGraph":
        newObjectGraph = buildObjectGraph()
        proj.objectGraph = newObjectGraph 
        proj.save_object_graph()
    elif command in ["createMain", "createSnakefile"]:
        proj.write_main_snakefile()
    elif command == "createSymbolicLinks":
        proj.create_symbolic_links()
    elif command in ["prepare","prepareTest"]:
        newObjectGraph = buildObjectGraph()
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
            proj.write_main_snakefile()
            proj.create_symbolic_links()
    elif command in ["prepareObjects"]:
        proj.create_symbolic_links()
    elif command == "run":
        sargs = ['snakemake',
                '-s', proj.get_pipeline_directory() + '/Snakefile', 
                '-d', proj.directory]
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

        default_remote_provider = get_arg_value(sargs,'--default-remote-provider')
        default_remote_prefix = get_arg_value(sargs,'--default-remote-prefix')
        if default_remote_provider and default_remote_prefix:
            upload_project_files_to_remote(default_remote_provider,
                                           default_remote_prefix)

        if ("--kubernetes" in sargs 
            or "--google-lifesciences" in sargs 
            or "--use-singularity" in sargs):
            if default_remote_provider and default_remote_prefix:
                os.environ['SO_REMOTE'] = f"{default_remote_provider}:{default_remote_prefix}"
            else:
                os.environ['SO_REMOTE'] = 'singularity'
            #sargs += ["--envvars SO_REMOTE "]
        os.execvp('snakemake',sargs)
    elif command == "submit":
        from snakeobjects.Project import ProjectException
        sargs = []
        if "default_snakemake_args" in proj.parameters:
            sargs += proj.parameters["default_snakemake_args"].split()
        else:
            raise ProjectException("No profile specified")
        sargs += args[1:]
        profile=sargs[sargs.index('--profile')+1]
        if not os.path.exists(profile): 
            raise ProjectException("Profile not found %s" % profile)
        if not os.path.exists(profile+"/config.yaml"): 
            raise ProjectException("No config.yaml in %s" % profile)
        pr_config = load_yaml(profile+"/config.yaml")
        if not "cluster" in pr_config: 
            ProjectException("cluster in not specified in %s" % profile+"/config.yaml")
        cmd=pr_config["cluster"]
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
        with open(proj.directory+'/objects/.snakeobjects/jobscript.sh', 'a') as js:
            for k,v in pr_config.items():
                if not k in 'jobname jobscript cluster cluster-status'.split(' '):
                    js.write('--'+str(k)+' ' + str(v) + ' ')
            js.write(' '.join(args[1:]))

        os.system("%s/%s" % (profile,cmd)+ " $SO_PROJECT/objects/.snakeobjects/jobscript.sh")
        #os.execvp('python', [profile + "/" +cmd, "$SO_PROJECT/objects/.snakeobjects/jobscript.sh"])        
    elif command == "describe":
        print("Project parameters:")
        for k,v in proj.parameters.items():
            print(f"\t{k}: {v}")
        proj.objectGraph.print_stats()
    elif command == "graph":
        print(args, file=sys.stderr)
        graph.driver(proj.objectGraph, args)
    elif command == "cleanProject":
        import shutil
        
        if os.path.exists(proj.directory+'/.snakemake'):
            val = input('Delete .snakemake? (y/n):')
            if val == 'y':
                shutil.rmtree(os.path.abspath(proj.directory+'/.snakemake'))
        if os.path.exists(proj.directory+'/OG.json'):
            val = input('Delete OG.json? (y/n):')
            if val == 'y':
                os.remove(os.path.abspath(proj.directory+'/OG.json')) 
        val = input('Delete all object directories? (y/n):')
        if val == 'y':   
            for ot in proj.objectGraph.get_object_types():
                if os.path.exists(proj.directory+'/'+ot):
                    shutil.rmtree(os.path.abspath(proj.directory+'/'+ot))
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
    #print("BBBBBBB")
    cli(sys.argv[1:])
