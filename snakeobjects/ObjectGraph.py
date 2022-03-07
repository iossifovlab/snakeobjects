import sys
import copy
import re
import os 
import json,yaml
from collections import OrderedDict

class OGO:
    """
    The class represents the objects in an :py:class:`.ObjectGraph`.
    Each object has the following attributes:

    .. py:attribute:: oType
       :type: str

       the object type of the object

    .. py:attribute:: oId
       :type: str

       the object id of the object

    .. py:attribute:: params 
       :type: dict[str,str] 

       the key-value parameters of the object

    .. py:attribute:: deps 
       :type: list[OGO] 

       the list of the dependency objects (other objects this object depends on)

    """
    def __init__(ogo,oType,oId,params=None,deps=None):
        ogo.oType = oType 
        ogo.oId = oId 
        ogo.params = params if params else {} 
        ogo.deps = deps if deps else []

        # for dep in ogo.deps:
        #     if isinstance(dep,(OGO)):
        #        dep.parents.append(ogo) 

    def k(o):
        return (o.oType,o.oId)

    def deepDeps(o, dot=None, level=1, mode='equal'):
        dps = []
        queue = []
        seen = set()
        out = []
        for d in o.deps:
            queue.append((d,1))
        while queue:
            x,xl = queue.pop(0)

            if x.k() in seen:
                if (mode == 'equal' and xl == level):
                    dps.append(x)
                continue
            seen.add(x.k())
            if ((mode == 'lessOrEqual' and xl <= level) or \
               (mode == 'equal' and xl == level)) or x.oType == dot:
                dps.append(x)
            for y in x.deps:
                queue.append((y,xl+1))

        # this is done to keep the order of objects, but output only unique
        for x in dps:
            if not x in out:
                out.append(x)
        if dot:
            out = [dp for dp in out if dp.oType == dot]
        return out

class ObjectGraphException(Exception):
    pass

class ObjectGraph:
    """
    The class representing the directed acyclic graph representing the relationships among
    the snakeobjects' objects. Each object is defined by a set an **object type** and by **object id**. 
    """

    def __init__(self):
        self.objectsByKey = {}
        self.O = {}
        self.oOrder = {}
        self.tOrder = []

    def __getitem__(self,key):
        """
        Query objects in the object graph. 

        :param key: the object query
        :type key: oType or (oType,oOid)

        The graph (G) indexing has two forms depending of the **key** parameter:

        1. ``G[oType]`` returns the **list** of the objects of type ``oType``;
        2. ``G[oType,oId]`` returns the object with object id ``oId`` and object type type ``oType``.

        :rtype: list(OGO) if form 1 is used or OGO if form 2 is used
        """
        if type(key) == tuple  and len(key) == 2:
            ot,oi = key
            return self.O[ot][oi]
        else:
            return self.oOrder[key]

    
    def get_object_types(self):
        """The object types used in the ObjectGraph"""
        
        return self.tOrder
    
    def add(self,oType,oId,params=None,deps=None):
        """
        Adds a new object to the object graph.
    
        :param str oType: the object type of the new object.
        :param str oId: the object id of the new object.
        :param params: the parameters of the new object. If None, the new 
                       object is assigned with no parameters.
        :type params: dict[str,str]
        :param deps: list objects the new object depends on. If None, the new object 
                     depends on no other objects. The objects included in deps should 
                     already be part of the ``ObjectGraph``.
        :type deps: list[OGO]
   
        :returns: the newly created object.
        :rtype: OGO
        """
        o = OGO(oType,oId,params,deps)
        self._addObject(o)
        return o
 
    def _addObject(self, ob):
        ob.OG = self
        key = str(ob.oType)+":"+str(ob.oId)

        if key in self.objectsByKey:
            raise ObjectGraphException(f"The object of type {ob.oType} and key {ob.oId} is already added\n");
        self.objectsByKey[key] = ob;
        if ob.oType in self.O:
            self.oOrder[ob.oType].append(ob)
            self.O[ob.oType][ob.oId] = ob
        else:
            self.O[ob.oType] = {}
            self.tOrder.append(ob.oType)
            self.O[ob.oType][ob.oId] = ob
            self.oOrder[ob.oType] = [ob]
                     
    def save(self,outFile):
        """
        Save the object graph in a json file ``outFile``.
        """
        f=open(outFile, "w") 

        sep = ""
        f.write("{\n")            
        for ot in self.tOrder:
            for o in self.oOrder[ot]:
                f.write(sep + "\t\""+str(o.oId)+"."+str(o.oType)+ "\": {\n")
                f.write("\t\t\"id\": \""+str(o.oId)+"\",\n")
                f.write("\t\t\"type\": \""+str(o.oType)+"\",\n")
                f.write("\t\t\"deps\":["  +  ",".join(["\""+d.oType+"/"+d.oId+"\"" for d in o.deps]) + "],\n" )
                f.write("\t\t\"params\": {\n")
                out = ""
                for p,v in sorted(o.params.items()):
                    w = "\""+str(v)+"\""
                    out += "\t\t\t\""+p+"\": "+w+",\n"
                f.write(out[:-2]+"\n")
                f.write("\t\t}\n")
                f.write("\t}")
                sep = ",\n\n"
        f.write("\n}\n")
        f.close()

    def print_stats(OG):
        print("Object types:")
        for tp,tpOs in sorted(OG.O.items()):
            print("\t", tp, ":", len(tpOs)) 
    

def load_object_graph(fname):
    with open(fname) as f:
        OD = json.load(f, object_pairs_hook=OrderedDict)
        OG = ObjectGraph()
        for k in OD:
            og = dict(OD[k])
            oi = og['id']
            ot = og['type']
            params = og['params'] if 'params' in og else []
            deps = og['deps'] if 'deps' in og else []
            dt = [d.split('/') for d in deps]
            OG.add(ot,oi, params, [OG[ot,oi] for ot,oi in dt])
        return OG

