import traceback 
from collections import OrderedDict

# _targetPrefix = "objLinks/"
_targetPrefix = ""
_config = None
_OG = None

def _find_object_type():
    fss = reversed(traceback.extract_stack())
    for fs in fss:
        if fs.filename.endswith(".snakefile"):
            ot = fs.filename.split("/")[-1][:-len(".snakefile")] 
            return ot
   
def set_object_graph(OG):
    global _OG
    _OG = OG 

def EFS(t):
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

def T(t): 
    return _targetPrefix + _find_object_type() + "/{oid}/"  + t

def TE(t): 
    return _targetPrefix + _find_object_type() + "/{{oid}}/"  + t

def P(p):
    ot = _find_object_type()
    return lambda wc: _OG[ot,wc.oid].params[p]

def PP(p):
    return _OG.params[p]

def all_obj_types():
    return _OG.tOrder

def DT(t, dot=None, level=1, mode='equal'): 
    ot = _find_object_type()
    def _DT(wc):
        dp = _OG[ot,wc.oid].deepDeps(dot,level,mode)
        return ["%s%s/%s/%s" % (_targetPrefix,d.type,d.name,t) for d in dp]
    return _DT 

def DP(p,dot=None, level=1, mode='equal'):
    ot = _find_object_type()
    def _DP(wc):
        dp = _OG[ot,wc.oid].deepDeps(dot,level,mode)
        return [d.params[p] for d in dp]
    return _DP
    
def all_obj_dirs(otype=None):
    if otype:
        return [_targetPrefix + otype + "/" + o.name for o in _OG[otype]]
    return [_targetPrefix + otype + "/" + o.name for otype in all_obj_types() for o in _OG[otype]]
