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

        .. py:attribute:: subprojects
           :type: dict[str,str] 

           a key value dictionary for subprojects names


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

        self.parent_projects = {}
        
        self._run_parameter_interpolation()

        if os.path.isfile(self._objectGraphFileName):
            self.objectGraph = load_object_graph(self._objectGraphFileName)
        else:
            self.objectGraph = ObjectGraph()

    def interpolate(self,O,oo=None):
        ptn = re.compile(r'(\[(\w+)\:([\w\/]+)(\:(\w+))?\])(\w*)')
        
        
        if type(O) == int or type(O) == float:
            return O
        elif type(O) == list:
            return [self.interpolate(x) for x in O]
        elif type(O) == dict:
            return {k:self.interpolate(v) for k,v in O.items()}
        elif type(O) == str:
            for s in ptn.findall(O):
                iType = s[1]
                name = s[2]
                #print("AAA",O,s)
                if iType == 'P':
                    if not oo:
                        raise ProjectException("object interpolation requires an object")
                    if name not in oo.params:
                        raise ProjectException("parameter %s is not defined for the object %s" % (name, oo.k()))
                    O = O.replace(s[0],oo.params[name])
                elif iType == 'E':
                    if not name in os.environ:
                        raise ProjectException("Varianble %s is not defined" % name)
                    O = O.replace(s[0],os.environ[name])
                elif iType == 'PP':
                    if not name in self.parameters:
                        raise ProjectException('Parameter %s is not defined' % name)
                    O = O.replace(s[0],self.parameters[name])
                elif iType == 'D':
                    if name == "project":
                        pv = self.directory
                    elif name == "pipeline":
                        pv = self.get_pipeline_directory()
                    else:
                        raise ProjectException('The project property %s is unknonw.' % name)
                    O = O.replace(s[0],pv)
                elif iType == 'NP':
                    if 'so_parent_projects' in self.parameters:
                        so_parent_projects = self.parameters['so_parent_projects']
                    else:
                        raise ProjectException('NP present, but there is no so_parents parameter')
                    if not '/' in name:
                        if name in so_parent_projects:
                            dir = so_parent_projects[name]
                            if not os.path.exists(dir):
                                raise ProjectException('The path to parent does not exist')
                            self.parent_projects[name] = Project(dir)
                            O = O.replace(s[0], str(self.parent_projects[name].parameters[s[4]]))
                        else:
                            raise ProjectException(f'Project {name} is not in so_parent_projects')
                    else:
                        cs = name.split('/')
                        if not cs[0] in so_parent_projects:
                            raise ProjectException(f'Project {cs[0]} is not in so_parent_projects')
                        else:
                            dir = so_parent_projects[cs[0]]
                            if not os.path.exists(dir):
                                raise ProjectException(f'The path {dir} to parent does not exist')
                            self.parent_projects[cs[0]] = Project(dir)
                            P = self.parent_projects[cs[0]]
                        for p in cs[1:]:
                            P = P.parent_projects[p]
                        O = O.replace(s[0], str(P.parameters[s[4]]))
                else:
                    raise ProjectException('Interpolation type [%s: ...] is unknown; can be only E|P|PP|D|NP.' % iType)
            
            return O

        
    def _run_parameter_interpolation(self):
        for k,v in self.parameters.items():
            self.parameters[k] = self.interpolate(v)
            

    def get_pipeline_directory(self):
        if "so_pipeline" in self.parameters:
            ppd = self.parameters['so_pipeline']
        elif "SO_PIPELINE" in os.environ:
            ppd = os.environ['SO_PIPELINE']
        elif os.path.exists("workflow"):
            return os.path.abspath("workflow")
        else:
            ppd = self.directory
        return os.path.abspath(ppd)

    def get_path(self,x='python'):
        path = ''
        name = 'so_extra_' + x + '_path'
        if name in self.parameters:
            v = self.parameters[name]
            path += ':'.join(v) if type(v) == list else v 
        for p in self.parent_projects.values():
            v = p.get_path(x)
            path += ':'+v if len(path) > 0 else v
        return path
        
    def save_object_graph(self):
        self.objectGraph.save(self._objectGraphFileName)

    def prepare_objects(self):
        self.write_main_snakefile()
        self.create_object_directories()
    
    def ensure_object_type_snakefile_exists(self,ot):
        sfile = self.get_pipeline_directory() + "/" + ot + ".snakefile"
        if not os.path.exists(sfile):
            print(f'WARNING: creating dummy snakefile {sfile.split("/")[-1]}')
            with open(sfile, 'w') as f:
                f.write(f'add_targets()\n')
        return sfile
        
    def write_main_snakefile(self):
        from glob import glob
        mf=self.get_pipeline_directory() + "/Snakefile"
        from snakeobjects import __version__
        o_types = set(self.objectGraph.get_object_types())
        o_list = set([x.split('/')[-1][:-10] for x in
                      glob(self.get_pipeline_directory() + '/*.snakefile')])
        all_o_types = sorted(set.union(o_types,o_list))

        if os.path.exists(mf):
            with open(mf) as f:
                old_version = f.readline().strip('\n\r').split(' ')[1]
                old_o_types = f.readline().strip('\n\r').split(' ')[2]
            if old_version == __version__ and old_o_types == ','.join(all_o_types):
                return 
        header = importlib_resources.read_text(__package__,'header.snakefile')
        with open(mf, 'w') as f:
            
            f.write(f'#snakeobjects {__version__}\n')
            f.write(f'#Snakeobjects types: {all_o_types}\n')
            f.write(header)

            for ot in all_o_types:
               
                sfile = self.ensure_object_type_snakefile_exists(ot) 
                sfile = ot + ".snakefile"
                f.write(f'include: "{sfile}"\n')

                f.write(f'rule so_all_{ot}:\n')
                f.write(f'  input:\n')
                f.write(f'    expand("{{of}}", of=project.get_all_object_flags("{ot}"))\n\n')

                f.write(f'rule so_{ot}_obj:\n')
                f.write(f'  input: get_targets("{ot}")\n')
                f.write(f'  output: touch("{ot}/{{oid}}/obj.flag")\n\n') 

    def get_object_flag(self,o):
        return f'{o.oType}/{o.oId}/obj.flag'

    def get_object_directory(self,o):
        return f'{self.directory}/{o.oType}/{o.oId}'

    def create_symbolic_links(self):
        for tp in sorted(self.objectGraph.get_object_types()):
            for o in self.objectGraph[tp]:
                oDir = self.get_object_directory(o)

                for k,v in list(o.params.items()):
                    if not k.startswith("symlink."):
                        continue
                    if not os.path.exists(oDir):
                        os.makedirs(oDir)
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
