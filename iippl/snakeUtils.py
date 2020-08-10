import traceback 

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
   
def set_config(config):
    global _config
    _config = config



def EFS(t):
    a = _targetPrefix + _find_object_type() + "/{oid}/log/" + t 
    r = {
        'E': a + '-err.txt',
        'O': a + '-out.txt',
        'T': a + '-time.txt'
    }
    print("EFS with", t, "returns", r)
    return r
 
def T(t): 
    return _targetPrefix + _find_object_type() + "/{oid}/"  + t

def DT(wc,t,dot=None): 
    # print("AAAA: DT called with =", wc, " and t=", t)
    ok = "%s.%s" % (wc.oid,_find_object_type())
    dp = _config[ok]["deps_local"]
    if dot:
        dp = [d for d in dp if d.startswith(dot)]
    r = ["%s%s/%s" % (_targetPrefix,d,t) for d in dp]
    # r = "objLinks/base/o/" + t
    # print("    : returning: ", r)
    return r
    

def all_obj_types():
    return {atts['type'] for atts in _config.values() if 'type' in atts}

def all_obj_ids(otype):
    # print('ALL TYPES:',all_obj_types())
    return sorted([x[:-(len(otype)+1)] for x in _config if x.endswith("." + otype)])

def all_obj_dirs(otype=None):
    if otype:
        return [_targetPrefix + otype + "/" + oid 
                    for oid in all_obj_ids(otype)]
    return [_targetPrefix + otype + "/" + oid 
                    for otype in all_obj_types() 
                            for oid in all_obj_ids(otype)]
