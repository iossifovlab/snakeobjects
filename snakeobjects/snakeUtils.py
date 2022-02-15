"""
"""

import traceback 
from collections import OrderedDict, defaultdict

# CLOUD related reorganization
_targetPrefix = ""
#_targetPrefix = ""
_OG = None
_project = None
_objectTypeTargets = defaultdict(list) 

def _find_object_type():
    fss = reversed(traceback.extract_stack())
    for fs in fss:
        if fs.filename.endswith(".snakefile"):
            ot = fs.filename.split("/")[-1][:-len(".snakefile")] 
            return ot

   
def set_project(project):
    global _OG,_project
    _project = project
    _OG = project.objectGraph


def add_targets(*targets):
    global _objectTypeTargets
    _objectTypeTargets[_find_object_type()] += targets


def get_targets(ot):
    global _objectTypeTargets, _OG

    def _GT(wc):
        return ["%s%s/%s/%s" % (_targetPrefix,d.oType,d.oId,"obj.flag")
                for d in _OG[ot,wc.oid].deps]
    return [_GT] + ["%s%s/{oid}/%s" % (_targetPrefix,ot,t)
                    for t in _objectTypeTargets[ot]] 


def PP(p,parent_project_id=None):
    """Value of the project parameter ``p`` in current project, if parent_project_id is not specified, otherwise the value of parameter ``p`` in ``parent_project_id``. If parameter ``p`` is not in the project parameters, then the value of this parameter in the first recusively searched parameters of parent projects. """
    return _project.get_parameter(p,parent_project_id)


def T(t): 
    """The file name for target ``t`` of the current object."""
    return _targetPrefix + _find_object_type() + "/{oid}/"  + t


def TE(t): 
    """
    The file name pattern for target ``t`` of the current object suitable as an
    argument to the expand function (a.k.a. the curly brackets are doubled).
    """
    return _targetPrefix + _find_object_type() + "/{{oid}}/"  + t


def EF(s): 
    """
    Argument ``s`` is a string that may include projects parameters as well as 
    object's parameters as they are definded with PP and P directives in the 
    so_project.yaml file.
    """
    ot = _find_object_type()
    def _DT(wc):
        o = _OG[ot,wc.oid]
        return _project.interpolate(s,o)
    return _DT 


def DT(t, dot=None, level=1, mode='equal'): 
    """
    List of files for the targets named ``t`` of the dependency objects (the objects the 
    current object depends on).

    :param str t: the name of the target in the dependency objects.

    :param str dot: the object type of the dependency objects to be used. If None, 
        dependency objects of all object types will be used.

    :param int level: how many levels down the dependency graph to explore. 
        The default value of 1 indicates that only the objects that the current object
        depends on directly will be used. 

    :param mode: if set to ``'equal'``, only objects that are exactly **level** steps away 
        from the current object will be used. If set to ``'lessOrEqual'``, objects of **level**
        steps or less will be used.
    :type mode: ``'equal'`` or ``'lessOrEqual'``
    
    """
    ot = _find_object_type()
    def _DT(wc):
        dp = _OG[ot,wc.oid].deepDeps(dot,level,mode)
        return ["%s%s/%s/%s" % (_targetPrefix,d.oType,d.oId,t) for d in dp]
    return _DT 


def P(p):
    """Value of the parameter ``p`` of the current object."""
    ot = _find_object_type()
    return lambda wc: _OG[ot,wc.oid].params[p]

def DP(p,dot=None, level=1, mode='equal'):
    """
    List of values of the parameter ``p`` of the dependency objects (the objects the 
    current object depends on).

    :param str p: the name of the parameter of the dependency objects.

    :param str dot: the object type of the dependency objects to be used. If None, 
        dependency objects of all object types will be used.

    :param int level: how many levels down the dependency graph to explore. 
        The default value of 1 indicates that only the objects that the current object
        depends on directly will be used. 

    :param mode: if set to ``'equal'``, only objects that are exactly **level** steps away 
        from the current object will be used. If set to ``'lessOrEqual'``, objects of **level**
        steps or less will be used.
    :type mode: ``'equal'`` or ``'lessOrEqual'``
    
    """
    ot = _find_object_type()
    def _DP(wc):
        dp = _OG[ot,wc.oid].deepDeps(dot,level,mode)
        return [d.params[p] for d in dp]
    return _DP
    

def LFS(t):
    """
    A set of three log files, log.O, log.E, and log.T, named after target ``t`` to be used in a rule. 

    This function is expected to be used in the ``log`` section of the snakemake rules 
    like ``log: **LFS('bamFile'))``. After that, the rule can use one or more of ``{log.O}``, ``{log.E}``, and 
    ``{log.T}`` for storing standard output, standard error, or timing information respectively. 
    """
    a = _targetPrefix + _find_object_type() + "/{oid}/log/" + t 
    r = {
        'E': a + '-err.txt',
        'O': a + '-out.txt',
        'T': a + '-time.txt'
    }
    return r


def B(t):
    a = _targetPrefix + _find_object_type() + "/{oid}/log/" + t 
    return a + '-bmk.txt'


