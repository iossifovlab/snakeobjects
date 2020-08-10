import traceback 

_targetPrefix = "objLinks/"
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
 
def T(t): 
    return _targetPrefix + _find_object_type() + "/{oid}/"  + t

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
