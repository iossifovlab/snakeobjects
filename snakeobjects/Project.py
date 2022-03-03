import os
import sys
import yaml
import re
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

class Pipeline:
    """
        Each objects of the Pipeline class represents a snakeobject pipeline. 
        
        A snakeobject pipeline attributes are:

        .. py:attribute:: snakefile_directory
           :type: str

           the snakefiles directory

        .. py:attribute:: environment_variables
           :type: dict[str,str] 

           a key value dictionary for environment variables

        .. py:attribute:: descrition
           :type: str

           a description of the pipeline

        Pipeline useful functions are:

        .. py:attribute:: get_snakefile_directory(object_instance)

           returns path to the Snakefile directory

        .. py:attribute:: get_environment_variables(object_instance)

           return returns dictionary with keys PATH, PYTHONPATH, etc.

        .. py:attribute:: get_description(object_instance )

           returns description of the pipeline
    """
        
    def __init__(self):
        self.snakefile_directory = None
        self.environment_variables = None    
        self.descriptions = "Pipeline contains location of snakefiles for the project"
        self.project = None
        
    def set_environment_variables(self):
        pass

    def get_snakefile_directory(self):
        return self.snakefile_directory
    
    def get_environment_variables(self):
        return self.environment_variables

    def get_description(self):
        return self.description

    
class DirectoryPipeline(Pipeline):
    """
    Subclass of the Pipeline class defined by directory path.
    """
    def __init__(self, dir, proj):
        super().__init__()
        self.snakefile_directory = dir
        self.project = proj
        self.set_environment_variables()
        
    def set_environment_variables(self):
        so_path = self.snakefile_directory
        path = [so_path]
        pythonpath = []
        perl5lib = []

        for x in 'bin python R perl'.split(' '):
            name = so_path + '/' + x
            if os.path.exists(name):
                path.append(name)
                if x == 'python':
                    pythonpath.append(name)
                if x == 'perl':
                    perl5lib.append(name)
            xname = 'so_extra_' + x + '_path'
            if xname in self.project.parameters:
                v = self.project.parameters[xname]
                path.append(v)
                if x == 'python':
                    pythonpath.append(v)

        for p in self.project.parent_projects.values():
            v = p.get_paths()

            if 'PATH' in v:
                path.append(v['PATH'])
            if 'PYTHONPATH' in v:
                pythonpath.append(v['PYTHONPATH'])
            if 'PERL5LIB' in v:
                perl5lib.append(v['PERL5LIB'])
                
        self.environment_variables = {'PATH': ':'.join(path),
                'PYTHONPATH': ':'.join(pythonpath),
                'PERL5LIB': ':'.join(perl5lib)
                }        
        
class PackagePipeline(Pipeline):
    """
    Subclass of Pipeline class defined by python package.
    """
    def __init__(self, config_string, proj):
        super().__init__()
        self.project = proj

        if not config_string.startswith('package:'):
            print(f"no package string {config_string}, should start with 'package:'")
            exit(1)
        else:
            import pkg_resources
            pkg_name, data_name = config_string.split(':')[1:]
            self.snakefile_directory=pkg_resources.resource_filename(pkg_name, data_name)

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

        .. py:attribute:: parent_projects
           :type: dict[str,str] 

           a key value dictionary for parent project names


        .. py:attribute:: objectGraph 
           :type: ObjectGraph

           the objectGraph for the project

        Project useful functions are:

        .. py:attribute:: get_object_directory(object_instance)

           returns directory of object

        .. py:attribute:: get_objects_path(type_name, object_name)

           return path to the objects of type `type_name`, and object name `object_name`

        .. py:attribute:: get_pipeline_directory( )

           returns path to the project pipeline directory

        .. py:attribute:: get_parameter(name, parent_project_id=None)

           If `parent_project_id` is given returns the value of parameter in that project, otherwise, if `name` is in project parameters, returns its value. Finally, if `name` is not in the project parameters, returns the value of named parameters in the first project recursively folowing parent projects. 

    """

    def __init__(self,directory=None):
        self.directory = directory if directory else find_so_project_directory()
        self._porojectFileName = self.directory + "/so_project.yaml"

        self._objectGraphFileName = self.directory + "/OG.json"

        if os.path.isfile(self._porojectFileName):
            self.parameters = load_yaml(self._porojectFileName)
            self.parameters['directory'] = self.directory
        else:
            self.parameters = {}

        self.parent_projects = {}
        
        self._run_parameter_interpolation()
        
        if 'so_parent_projects' in self.parameters and not self.parent_projects:
            for p,d in self.parameters['so_parent_projects'].items():
                self.parent_projects[p] = Project(d)

        if os.path.isfile(self._objectGraphFileName):
            self.objectGraph = load_object_graph(self._objectGraphFileName)
        else:
            self.objectGraph = ObjectGraph()

        self.set_pipeline_directory()

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
                prnt = s[4]
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
                elif iType == 'D':
                    if name == "project":
                        pv = self.directory
                    elif name == "pipeline":
                        pv = self.get_pipeline_directory()
                    else:
                        raise ProjectException('The project property %s is unknonw.' % name)
                    O = O.replace(s[0],pv)
                elif iType == 'PP':
                    if ('so_parent_projects' in self.parameters
                        and not self.parent_projects):
                        for k,v in self.parameters['so_parent_projects'].items():
                            self.parent_projects[k] = Project(v)                    
                    if not prnt:
                        if not name in self.parameters:
                            raise ProjectException('Parameter %s is not defined' % name)
                        O = O.replace(s[0],self.parameters[name])
                    elif not '/' in prnt:
                        if not '/' in prnt:
                            if prnt in self.parent_projects:
                                dir = self.parent_projects[prnt].directory
                                if not os.path.exists(dir):
                                    raise ProjectException('The path to parent does not exist')
                                self.parent_projects[prnt] = Project(dir)
                                O = O.replace(s[0], str(self.parent_projects[prnt].parameters[name]))
                            else:
                                raise ProjectException(f'Project {prnt} is not in so_parent_projects')
                    elif '/' in prnt:
                            cs = prnt.split('/')
                            if not cs[0] in self.parent_projects:
                                raise ProjectException(f'Project {cs[0]} is not in so_parent_projects')
                            else:
                                dir = self.parent_projects[cs[0]]
                                if not os.path.exists(dir):
                                    raise ProjectException(f'The path {dir} to parent does not exist')
                                self.parent_projects[cs[0]] = Project(dir)
                                P = self.parent_projects[cs[0]]
                            for p in cs[1:]:
                                P = P.parent_projects[p]
                            O = O.replace(s[0], str(P.parameters[name]))
                    else:
                        raise ProjectException('Interpolation type [%s: ...] is unknown; can be only E|P|PP|D.' % iType)
            
            return O

    def _run_parameter_interpolation(self):
        for k,v in sorted(self.parameters.items(), key=lambda x: not '[D:' in x):
            self.parameters[k] = self.interpolate(v)
        
    def get_parent_project(self, parent_project_id):
        parent_project_ids = parent_project_id.split("/")
        proj = self
        for id in parent_project_ids:
            proj = proj.parent_projects[id]
        return proj

    def get_parameter(self, name, parent_project_id=None):
        if parent_project_id:
            proj = self.get_parent_project(parent_project_id)
            return proj.parameters[name] 
        if name in self.parameters:
            return self.parameters[name]
        for proj in self.parent_projects.values():
            v = proj.get_parameter(name)
            if v:
                 return v 
        return None

    def set_pipeline_directory(self):
        if "so_pipeline" in self.parameters:
            ppd = self.parameters['so_pipeline']
            if ppd.startswith('directory:'):
                ppd = os.path.abspath(ppd.split(':')[1])
                self.pipeline = DirectoryPipeline(ppd, self)
            elif ppd.startswith('package'):
                self.pipeline = PackagePipeline(ppd, self)
            else:
                self.pipeline = DirectoryPipeline(os.path.abspath(ppd), self)
        elif "SO_PIPELINE" in os.environ:
            self.pipeline = DirectoryPipeline(os.path.abs.path(os.environ['SO_PIPELINE']))
        elif os.path.exists("workflow"):
            self.pipeline = DirectoryPipeline(os.path.abspath("workflow"), self)
        else:
            self.pipeline = DirectoryPipeline(os.path.abspath(self.directory))
            
    def get_pipeline_directory(self):
        
        if "so_pipeline" in self.parameters:
            ppd = self.parameters['so_pipeline']
            if ppd.startswith('directory:'):
                ppd = os.path.abspath(ppd.split(':')[1])
            elif ppd.startswith('package'):
                pclass = PackagePipeline(ppd, self)
                ppd = pclass.get_snakefile_directory()
            else:
                ppd = os.path.abspath(ppd)
        elif "SO_PIPELINE" in os.environ:
            ppd = os.path.abs.path(os.environ['SO_PIPELINE'])
        elif os.path.exists("workflow"):
            ppd = os.path.abspath("workflow")
        else:
            ppd = os.path.abspath(self.directory)
        if os.path.exists(ppd):
            return ppd
        else:
            print(f"so_pipeline file {ppd} does not exists")
            exit(1)

    def get_paths(self):
        so_path = self.get_pipeline_directory()
        path = [so_path]
        pythonpath = []
        perl5lib = []

        for x in 'bin python R perl'.split(' '):
            name = so_path + '/' + x
            if os.path.exists(name):
                path.append(name)
                if x == 'python':
                    pythonpath.append(name)
                if x == 'perl':
                    perl5lib.append(name)
            xname = 'so_extra_' + x + '_path'
            if xname in self.parameters:
                v = self.parameters[xname]
                path.append(v)
                if x == 'python':
                    pythonpath.append(v)

        for p in self.parent_projects.values():
            v = p.get_paths()

            if 'PATH' in v:
                path.append(v['PATH'])
            if 'PYTHONPATH' in v:
                pythonpath.append(v['PYTHONPATH'])
            if 'PERL5LIB' in v:
                perl5lib.append(v['PERL5LIB'])
                
        return {'PATH': ':'.join(path),
                'PYTHONPATH': ':'.join(pythonpath),
                'PERL5LIB': ':'.join(perl5lib)
                }
        
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

    def get_object_path(self,oType,oId):
        return Path(self.get_object_directory(self.objectGraph[oType,oId]))

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
