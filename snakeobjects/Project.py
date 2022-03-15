from __future__ import annotations
import importlib
import os
import re
from typing import Dict, List
from pathlib import Path
from abc import ABC
import importlib.resources as importlib_resources
from importlib.util import spec_from_file_location, module_from_spec
import yaml
from collections import defaultdict

from snakeobjects.ObjectGraph import ObjectGraph, load_object_graph


def load_yaml(file_name):

    CF = open(file_name, 'r')
    config = yaml.safe_load(CF)
    CF.close()

    return config if config else {}


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


class Pipeline(ABC):
    """
        Each objects of the Pipeline class represents a snakeobject pipeline.
    """

    def get_snakefile_directory(self) -> Path:
        raise Exception("Should not be called!")

    def get_main_snakefile_path(self) -> Path:
        raise Exception("Should not be called!")

    def get_environment_variables(self) -> Dict[str, List[str]]:
        raise Exception("Should not be called!")

    def get_definition(self) -> str:
        '''
        Definitions is a string of the form:
             "directory:<dir>" OR
             "package:<package>" OR
             "<dir>" which is equivallent to "directory:<dir>"
        '''
        raise Exception("Should not be called!")

    def build_object_graph(self, project: Project, args: List[str]) \
                -> ObjectGraph:
        project = project
        args = args
        raise Exception("not yet")


class DirectoryPipeline(Pipeline):
    """
    Subclass of the Pipeline class defined by directory path.
    """

    def __init__(self, definition_dir: str):
        self.definition_dir = definition_dir
        self.snakefile_directory = Path(definition_dir).resolve(True)

    def get_snakefile_directory(self) -> Path:
        return self.snakefile_directory

    def get_main_snakefile_path(self) -> Path:
        return self.snakefile_directory / "Snakefile"

    DIR_TO_ENV = {
        "bin": ["PATH"],
        "python": ["PATH", "PYTHONPATH"],
        "R": ["PATH"],
        "perl": ["PATH", "PERL5LIB"]
    }

    def get_environment_variables(self) -> Dict[str, List[str]]:
        vars: Dict[str, List[str]] = {}
        vars["PATH"] = [str(self.snakefile_directory)]
        for sub_dir, env_vars in DirectoryPipeline.DIR_TO_ENV.items():
            sub_dir_path = self.snakefile_directory / sub_dir
            if os.path.exists(sub_dir_path):
                for var in env_vars:
                    try:
                        vars[var].append(str(sub_dir_path))
                    except KeyError:
                        vars[var] = [str(sub_dir_path)]
        return vars

    def get_definition(self) -> str:
        return self.definition_dir
        # return f"directory:{self.definition_dir}"

    def build_object_graph(self, project: Project, args: List[str]) \
            -> ObjectGraph:
        bldObjGraphPy = self.get_snakefile_directory() / "build_object_graph.py"
        if not os.path.isfile(bldObjGraphPy):
            raise Exception(f'ERROR: There is no {bldObjGraphPy}')

        spec = spec_from_file_location("build_object_graph", bldObjGraphPy)
        if spec is None:
            raise Exception("Could build specs for build_object_graph")
        if spec.loader is None:
            raise Exception("Loader is empty??")
        foo = module_from_spec(spec)
        spec.loader.exec_module(foo)
        newObjectGraph = ObjectGraph()
        foo.run(project, newObjectGraph, *args)
        return newObjectGraph


class PackagePipeline(Pipeline):
    """
    Subclass of Pipeline class defined by python package.
    """

    def __init__(self, definition_package: str):
        self.definition_package = definition_package
        self.pipeline_package = importlib.import_module(self.definition_package)

        if self.pipeline_package.__file__ is None:
            raise Exception(f"{definition_package} should be a python package 1")
        if not self.pipeline_package.__file__.endswith('__init__.py'):
            raise Exception(f"{definition_package} should be a python package 2")

        self.snake_file_dir = Path(self.pipeline_package.__file__).parent

    def get_snakefile_directory(self) -> Path:
        return self.snake_file_dir

    def get_main_snakefile_path(self) -> Path:
        return self.snake_file_dir / "Snakefile"

    def get_environment_variables(self) -> Dict[str, List[str]]:
        env = {'python':['PATH', 'PYTHONPATH'], 'bin':'PATH', 'perl':['PATH','PERL5LIB']}
        vars = defaultdict(list)
        for key,value in env.items():
            pa = self.snake_file_dir / key
            if pa.resolve().exists():
                if type(value) == list:
                    for v in value:
                        vars[v].append(str(pa))
                else:
                    vars[value].append(str(pa))
        return vars

    def get_definition(self) -> str:
        return f"pacakage:{self.definition_package}"

    def build_object_graph(self, project: Project, args: List[str]) \
            -> ObjectGraph:
        new_object_graph = ObjectGraph()
        self.pipeline_package.build_object_graph(
            project, new_object_graph, args)
        return new_object_graph


def build_pipeline(pipeline_definition: str) -> Pipeline:
    if pipeline_definition.startswith('directory:'):
        dir_path = pipeline_definition.split(':')[1]
        return DirectoryPipeline(dir_path)

    if pipeline_definition.startswith('package'):
        package_name = pipeline_definition.split(':')[1]
        return PackagePipeline(package_name)

    return DirectoryPipeline(pipeline_definition)


class Project:
    """
        Each objects of the Project class represents a snakeobject project.

        The project directory is given as the ``directory`` parameter. If 
        ``directory`` is none,
        the :py:func:`find_so_project_directory` function is used to determine 
        the project directory.

        A snakeobject project attributes are:

        .. py:attribute:: directory
           :type: str

           the project directory

        .. py:attribute:: parameters
           :type: dict[str,str]

           a key value dictionary for global project level parameters

        .. py:attribute:: parent_projects
           :type: dict[str,Project]

           a key value dictionary for parent project names


        .. py:attribute:: objectGraph
           :type: ObjectGraph

           the objectGraph for the project

        Project useful functions are:

        .. py:attribute:: get_object_directory(object_instance)

           returns directory of object

        .. py:attribute:: get_object_path(type_name, object_name)

           return path to the object of type `type_name`, and object name `object_name`

        .. py:attribute:: get_pipeline_directory( )

           returns path to the project pipeline directory

        .. py:attribute:: get_parameter(name, parent_project_id=None)

           If `parent_project_id` is given returns the value of parameter in
           that project, otherwise, if `name` is in project parameters, returns
           its value. Finally, if `name` is not in the project parameters,
           returns the value of named parameters in the first project
           recursively folowing parent projects.

    """

    def __init__(self, directory=None):
        self.directory = directory if directory else find_so_project_directory()
        self._porojectFileName = self.directory + "/so_project.yaml"

        self._objectGraphFileName = self.directory + "/OG.json"

        if os.path.isfile(self._porojectFileName):
            self.parameters = load_yaml(self._porojectFileName)
            self.parameters['directory'] = self.directory
        else:
            self.parameters = {}

        self.parent_projects: Dict[str, Project] = {}

        self._run_parameter_interpolation()

        if 'so_parent_projects' in self.parameters and not self.parent_projects:
            for p, d in self.parameters['so_parent_projects'].items():
                self.parent_projects[p] = Project(d)

        if os.path.isfile(self._objectGraphFileName):
            self.objectGraph = load_object_graph(self._objectGraphFileName)
        else:
            self.objectGraph = ObjectGraph()

        self.set_pipeline()

    def interpolate(self, obj, oo=None):
        ptn = re.compile(r'(\[(\w+)\:([\w\/]+)(\:(\w+))?\])(\w*)')

        if type(obj) == int or type(obj) == float:
            return obj
        elif type(obj) == list:
            return [self.interpolate(x) for x in obj]
        elif type(obj) == dict:
            return {k: self.interpolate(v) for k, v in obj.items()}
        elif type(obj) == str:
            for s in ptn.findall(obj):
                iType = s[1]
                name = s[2]
                prnt = s[4]
                # print("AAA",O,s)
                if iType == 'P':
                    if not oo:
                        raise ProjectException("object interpolation requires an object")
                    if name not in oo.params:
                        raise ProjectException(
                            "parameter %s is not defined for the object %s" % (name, oo.k()))
                    obj = obj.replace(s[0], oo.params[name])
                elif iType == 'E':
                    if name not in os.environ:
                        raise ProjectException("Varianble %s is not defined" % name)
                    obj = obj.replace(s[0], os.environ[name])
                elif iType == 'D':
                    if name == "project":
                        pv = self.directory
                    elif name == "pipeline":
                        pv = str(self.get_pipeline_directory())
                    else:
                        raise ProjectException('The project property %s is unknonw.' % name)
                    obj = obj.replace(s[0], pv)
                elif iType == 'PP':
                    if ('so_parent_projects' in self.parameters and
                            not self.parent_projects):
                        for k, v in self.parameters['so_parent_projects'].items():
                            self.parent_projects[k] = Project(v)
                    if not prnt:
                        if name not in self.parameters:
                            raise ProjectException('Parameter %s is not defined' % name)
                        obj = obj.replace(s[0], self.parameters[name])
                    elif '/' not in prnt:
                        if '/' not in prnt:
                            if prnt in self.parent_projects:
                                dir = self.parent_projects[prnt].directory
                                if not os.path.exists(dir):
                                    raise ProjectException('The path to parent does not exist')
                                self.parent_projects[prnt] = Project(dir)
                                obj = obj.replace(
                                    s[0], str(self.parent_projects[prnt].parameters[name]))
                            else:
                                raise ProjectException(
                                    f'Project {prnt} is not in so_parent_projects')
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
                        obj = obj.replace(s[0], str(P.parameters[name]))
                    else:
                        raise ProjectException(
                            'Interpolation type [%s: ...] is unknown; can be only E|P|PP|D.' % iType)

            return obj

    def _run_parameter_interpolation(self):
        for k, v in sorted(self.parameters.items(), key=lambda x: not '[D:' in x):
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

    def set_pipeline(self):
        pipeline_definition = f"directory:{self.directory}"
        if "SO_PIPELINE" in os.environ:
            pipeline_definition = os.environ['SO_PIPELINE']
        if "so_pipeline" in self.parameters:
            pipeline_definition = self.parameters["so_pipeline"]
        self.pipeline = build_pipeline(pipeline_definition)

        '''
        if "so_pipeline" in self.parameters:
            ppd = self.parameters['so_pipeline']
            if ppd.startswith('directory:'):
                ppd = os.path.abspath(ppd.split(':')[1])
                self.pipeline = DirectoryPipeline(ppd)
            elif ppd.startswith('package'):
                self.pipeline = PackagePipeline(ppd.split(':')[1])
            else:
                self.pipeline = DirectoryPipeline(ppd)
        elif "SO_PIPELINE" in os.environ:
            self.pipeline = DirectoryPipeline(os.environ['SO_PIPELINE'])
        elif os.path.exists("workflow"):
            self.pipeline = DirectoryPipeline(os.path.abspath("workflow"))
        else:
            self.pipeline = DirectoryPipeline(os.path.abspath(self.directory))
        '''

    def get_pipeline_directory(self) -> Path:
        self.set_pipeline()
        return self.pipeline.get_snakefile_directory()

    def get_environment_variables(self) -> Dict[str, List[str]]:
        #vars: Dict[str, List[str]] = {}
        vars = defaultdict(list)
        def add_vars(extra_vars: Dict[str, List[str]]):
            for var, values in extra_vars.items():
                try:
                    vars[var] += values
                except KeyError:
                    vars[var] = values

        # so_environment_EEE_add
        add_env_re = re.compile("so_environment_(.*)_add")
        for param_key, param_value in self.parameters.items():
            match = add_env_re.match(param_key)
            if match:
                env_name = match.group(1)
                if type(param_value) == list:
                    for pv in param_value:
                        vars[env_name].append(pv)
                else:
                    vars[env_name].append(param_value)

        add_vars(self.pipeline.get_environment_variables())

        for pp in self.parent_projects.values():
            add_vars(pp.get_environment_variables())

        # so_environment_EEE_set
        set_env_re = re.compile("so_environment_(.*)_set")
        for param_key, param_value in self.parameters.items():
            match = set_env_re.match(param_key)
            if match:
                env_name = match.group(1)
                vars[env_name] = [param_value]
        return vars

    def save_object_graph(self):
        self.objectGraph.save(self._objectGraphFileName)

    def ensure_object_type_snakefile_exists(self, ot):
        sfile = self.get_pipeline_directory() / (ot + ".snakefile")
        if not os.path.exists(sfile):
            print(f'WARNING: creating dummy snakefile '
                  f'{str(sfile).split("/")[-1]}')
            with open(sfile, 'w') as f:
                f.write(f'add_targets()\n')
        return sfile

    def write_main_snakefile(self):
        from glob import glob
        mf = self.get_pipeline_directory() / "Snakefile"
        from snakeobjects import __version__
        o_types = set(self.objectGraph.get_object_types())
        o_list = set([x.split('/')[-1][:-10] for x in
                      glob(str(self.get_pipeline_directory() / '*.snakefile'))])
        all_o_types = sorted(set.union(o_types, o_list))

        if os.path.exists(mf):
            with open(mf) as f:
                old_version = f.readline().strip('\n\r').split(' ')[1]
                old_o_types = f.readline().strip('\n\r').split(' ')[2]
            if old_version == __version__ and old_o_types == ','.join(all_o_types):
                return
        header = importlib_resources.read_text(__package__, 'header.snakefile')
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

    def get_object_flag(self, o):
        return f'{o.oType}/{o.oId}/obj.flag'

    def get_object_directory(self, o):
        return f'{self.directory}/{o.oType}/{o.oId}'

    def get_object_path(self, oType, oId):
        return Path(self.get_object_directory(self.objectGraph[oType, oId]))

    def create_symbolic_links(self):
        for tp in sorted(self.objectGraph.get_object_types()):
            for o in self.objectGraph[tp]:
                oDir = self.get_object_directory(o)

                for k, v in list(o.params.items()):
                    if not k.startswith("symlink."):
                        continue
                    if not os.path.exists(oDir):
                        os.makedirs(oDir)
                    dst = oDir + "/" + k[8:]
                    src = v
                    os.system("ln -sf %s %s" % (src, dst))
                    os.system("touch -h -r %s %s" % (src, dst))

    def get_all_object_flags(self, oType=None):
        OG = self.objectGraph
        if oType:
            return [self.get_object_flag(o) for o in OG[oType]]
        return [self.get_object_flag(o) for oType in OG.get_object_types() for o in OG[oType]]

    def set_environment(self):
        print("UPDATING ENVIRONMENT:")
        print("export SO_PROJECT=", self.directory, sep="")
        print("export SO_PIPELINE=", self.pipeline.get_definition(), sep="")

        os.environ['SO_PROJECT'] = self.directory
        os.environ['SO_PIPELINE'] = str(self.pipeline.get_definition())

        vars = self.get_environment_variables()

        for var, values in vars.items():
            my_val = ":".join(values)
            if var in ['PATH', 'PYTHONPATH', 'PERL5LIB'] and var in os.environ:
                print(f"export {var}={my_val}:${var}")
                values.append(os.environ[var])
            else:
                print(f"export {var}={my_val}")
            full_val = ":".join(values)
            os.environ[var] = full_val

    def buildObjectGraph(self, args: List[str]):
        self.set_environment()
        self.objectGraph = self.pipeline.build_object_graph(self, args)


if __name__ == "__main__":
    print("hi")
    p = Project()
    print("Project directory", p.directory)
    print("Project parameters", p.parameters)
    print("Number of object types in object graph", len(p.objectGraph.O))
    print("Number of objects in object graph",
        sum([len(x) for x in p.objectGraph.O.values()]))
    print("Pipeline directory is", p.pipeline.get_definition())
