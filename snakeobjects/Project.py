import os,sys,yaml,re
from pathlib import Path
import importlib.resources as importlib_resources
from snakeobjects.ObjectGraph import ObjectGraph, load_object_graph

def load_yaml(file_name):

    CF = open(file_name, 'r')
    config = yaml.safe_load(CF)
    CF.close()

    return config  if config else {}


def find_so_project_directory():
    """     
    Determines the directory for the current snakeobjects project.
    
    Several approaches are attempted and the first suitable directory found 
    is returned. If everything fails, the current working directory is returned.

    1. If the SO_PROJECT environment variable exits, its values is returned.
    2. If the current working directory contains a file named so_project.yaml, 
       the current working directory is returned.
    3. The parents of the current working directory are examined recursively 
       and the first one containing so_project.yaml file is returned. 
    """
    if "SO_PROJECT" in os.environ:
        return os.environ["SO_PROJECT"]
    cw = Path(os.getcwd())
    for d in [cw] + list(cw.parents):
        if os.path.isfile(d / "so_project.yaml"):
            return str(d)
    return os.getcwd()

class ProjectException(Exception):
    pass

class Project:
    """
        Each objects of the Project class represents a snakeobject project. 
        
        The project directory is given as the ``directory`` parameter. If ``directory`` is none,
        the :py:func:`find_so_project_directory` function is used to determine the project directory.

        A snakeobject project attributes are:

        .. py:attribute:: directory 
           :type: str

           the project directory

        .. py:attribute:: parameters 
           :type: dict[str,str] 

           a key value dictionary for global project level parameters

        .. py:attribute:: objectGraph 
           :type: ObjectGraph

           the objectGraph for the project
    """

    def __init__(self,directory=None):
        self.directory = directory if directory else find_so_project_directory()
        self._porojectFileName = self.directory + "/so_project.yaml"

        self._objectGraphFileName = self.directory + "/OG.json"

        if os.path.isfile(self._porojectFileName):
            self.parameters = load_yaml(self._porojectFileName)
        else:
            self.parameters = {}

        self._run_parameter_interpolation()

        if os.path.isfile(self._objectGraphFileName):
            self.objectGraph = load_object_graph(self._objectGraphFileName)
        else:
            self.objectGraph = ObjectGraph()

    def _run_parameter_interpolation(self):

        def interpolate(O):
            ptn = re.compile(r'(\[[\w]\:(\w+)\])(\w*)')
            
            if type(O) == int or type(O) == float:
                return O
            elif type(O) == list:
                return [interpolate(x) for x in O]
            elif type(O) == dict:
                return {k:interpolate(v) for k,v in O.items()}
            elif type(O) == str:
                for s in ptn.findall(O):
                    letter = s[0][1]
                    name = s[1]
                    if letter =='E':
                        if not name in os.environ:
                            raise ProjectException("Varianble %s is not defined" % name)
                        O = O.replace(s[0],os.environ[name])
                    elif letter =='C':
                        if not name in self.parameters:
                            raise ProjectException('Parameter %s is not defined' % name)
                        O = O.replace(s[0],self.parameters[name])
                    elif letter =='P':
                        if name == "projectDir":
                            pv = self.directory
                        elif name == "pipelineDir":
                            pv = self.get_pipeline_directory()
                        else:
                            raise ProjectException('The project property %s is unknonw.' % name)
                        O = O.replace(s[0],pv)
                    else:
                        raise ProjectException('Letter [%s: ...] is unknown; can be only E|C|P.' % letter)
                return O
            
        for k,v in self.parameters.items():
            self.parameters[k] = interpolate(v)
            

    def get_pipeline_directory(self):
        if "so_pipeline" in self.parameters:
            ppd = self.parameters['so_pipeline']
            #if not os.path.isabs(ppd):
            #    ppd = self.directory + "/" + ppd
        elif "SO_PIPELINE" in os.environ:
            ppd = os.environ['SO_PIPELINE']
        else:
            ppd = self.directory
        #return os.path.abspath(ppd)
        return ppd


    '''
    def prepare(self, newObjectGraph, ARGV=None):

        if not ARGV: ARGV = sys.argv 
        cmd = 'prepareTest' if len(ARGV) < 2 else ARGV[1]

        if cmd == 'prepareTest':
            print("Current graph stats")
            print("+++++++++++++++++++")
            self.objectGraph.print_stats()
            print("\n")
            print("New graph stats")
            print("+++++++++++++++")
            newObjectGraph.print_stats()
        elif cmd == 'prepare':
            self.objectGraph = newObjectGraph 
            os.system("mkdir -p objects/.snakeobjects")
            self.objectGraph.save(self.directory + "/objects/.snakeobjects/OG.json")
            self.prepare_objects()
        else:
            print(f"unkown command {cmd}. The known commands are 'prepareTest' and 'prepare'")
    '''
    def ensure_snakeobject_private_directory(self):
        #sopd = self.directory + "/objects/.snakeobjects"
        #if not os.path.exists(sopd):
        #    os.makedirs(sopd)
        #return sopd
        return self.directory
    
    def save_object_graph(self):
        self.objectGraph.save(self.ensure_snakeobject_private_directory() + "/OG.json")

    def prepare_objects(self):
        self.write_main_snakefile()
        self.create_object_directories()

    def ensure_object_type_snakefile_exists(self,ot):
        #sfile = self.get_pipeline_directory() + "/" + ot + ".snakefile"
        sfile = ot + ".snakefile"
        if not os.path.exists(sfile): 
            with open(sfile, 'w') as f:
                f.write(f'add_targets()\n')
        return sfile
        
    def write_main_snakefile(self):
        mf=self.ensure_snakeobject_private_directory() + "/Snakefile"
        header = importlib_resources.read_text(__package__,'header.snakefile')
        with open(mf, 'w') as f:
            f.write(header)

            for ot in self.objectGraph.get_object_types():
               
                sfile = self.ensure_object_type_snakefile_exists(ot) 

                f.write(f'include: "{sfile}"\n')

                f.write(f'rule so_all_{ot}:\n')
                f.write(f'  input:\n')
                f.write(f'    expand("{{of}}", of=project.get_all_object_flags("{ot}"))\n\n')

                f.write(f'rule so_{ot}_obj:\n')
                f.write(f'  input: get_targets("{ot}")\n')
                f.write(f'  output: touch("objects/{ot}/{{oid}}/obj.flag")\n\n') 

    def get_object_flag(self,o):
        return f'objects/{o.oType}/{o.oId}/obj.flag'

    def get_object_directory(self,o):
        return f'{self.directory}/objects/{o.oType}/{o.oId}'

    def create_object_directories(self):
        for tp in sorted(self.objectGraph.get_object_types()):
            for o in self.objectGraph[tp]:
                oDir = self.get_object_directory(o)

                logDir = oDir + "/log"
                if not os.path.exists(logDir):
                    os.makedirs(logDir)

                for k,v in list(o.params.items()):
                    if not k.startswith("symlink."):
                        continue
                    dst = oDir + "/" + k[8:]
                    src = v
                    os.system("ln -sf %s %s" % (src,dst))
                    os.system("touch -h -r %s %s" % (src,dst))

    def get_all_object_flags(self,oType=None):
        OG = self.objectGraph
        if oType:
            return [self.get_object_flag(o) for o in OG[oType]]
        return [self.get_object_flag(o) for oType in OG.get_object_types() for o in OG[oType]]

if __name__ == "__main__":
    print("hi")
    p = Project()
    print("Project directory",p.directory)
    print("Project parameters",p.parameters)
    print("Number of object types in object graph",len(p.objectGraph.O))
    print("Number of objects in object graph",sum([len(x) for x in p.objectGraph.O.values()]))
    print("Pipeline directory is",p.get_pipeline_directory())
