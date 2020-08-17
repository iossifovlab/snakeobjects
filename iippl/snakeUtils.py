import traceback 
from collections import OrderedDict

# _targetPrefix = "objLinks/"
_targetPrefix = ""
_config = None

def _find_object_type():
    fss = reversed(traceback.extract_stack())
    for fs in fss:
        if fs.filename.endswith(".snakefile"):
            ot = fs.filename.split("/")[-1][:-len(".snakefile")] 
            # print("FOUND OT", ot)
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

 
def P(p):
    ot = _find_object_type()
    return lambda wc: return _OG[ot,wc.oid].param[p]

def T(t): 
    return _targetPrefix + _find_object_type() + "/{oid}/"  + t

def all_obj_types():
    return _OG.tOrder

def deps(dp, level=0, mode=False):
    if level == 0:
        return dp
    if mode:
        for d in dp:
            ob = ".".join(d.split("/")[::-1])
            dp += deps(_config[ob]["deps_local"] 
                       if "deps_local" in _config[ob] else [], 
                       level - 1, mode)
        return dp

def DT(t,dot=None, level=1, mode=False): 
    # print("AAAA: DT called with =", wc, " and t=", t)
    ot = _find_object_type()
    def _DT(wc):
        ok = "%s.%s" % (wc.oid,ot)
        dp = _config[ok]["deps_local"] if "deps_local" in _config[ok] else []
        dp = deps(dp, level-1, mode)
        dp = list(OrderedDict.fromkeys(dp))
        if dot:
            dp = [d for d in dp if d.startswith(dot)]
        r = ["%s%s/%s" % (_targetPrefix,d,t) for d in dp]
        # r = "objLinks/base/o/" + t
        # print("    : returning: ", r)
        return r
    return _DT 

def DP(p,dot=None):
    ot = _find_object_type()
    def _DP(wc):
        ok = "%s.%s" % (wc.oid,ot)
        dp = _config[ok]["deps_local"] if "deps_local" in _config[ok] else []
        if dot:
            dp = [d for d in dp if d.startswith(dot)]
        r = []
        for d in dp:
            dok = ".".join((reversed(d.split("/"))))
            # print("DDDDDD",d,dok)
            r.append(_config[dok]['params'][p])
        return r
    return _DP
    

def all_obj_dirs(otype=None):
    if otype:
        return [_targetPrefix + otype + "/" + oid for oid in _OG[otype]]
    return [_targetPrefix + otype + "/" + oid for otype in all_obj_types() for oid in _OG[otype]]
