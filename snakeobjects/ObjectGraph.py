#!/usr/bin/env python


import sys
import copy
import re
import os 
import json,yaml
from collections import OrderedDict

class OGO:
    def __init__(ogo,tp,name,params=None,deps=None):
        ogo.type = tp 
        ogo.name = name 
        ogo.params = params if params else {} 
        ogo.deps = deps if deps else []
        ogo.parents = []
        ogo.OG = None
        for dep in ogo.deps:
            if isinstance(dep,(OGO)):
                dep.parents.append(ogo) 

    def k(o):
        return (o.type,o.name)

    def deepDeps(o, dot=None, level=1, mode='equal'):
        dps = []
        queue = []
        seen = set()
        for d in o.deps:
            queue.append((d,1))
        while queue:
            x,xl = queue.pop(0)

            if x.k() in seen: continue
            seen.add(x.k())
            if (mode == 'lessOrEqual' and xl <= level) or \
               (mode == 'equal' and xl == level):
                dps.append(x)
            for y in x.deps:
                queue.append((y,xl+1)) 
        if dot:
            dps = [dp for dp in dps if dp.type == dot]    
        return dps

class ObjectGraphException(Exception):
    pass

class ObjectGraph:

    def __init__(self):
        self.objectsByKey = {}
        self.O = {}
        self.oOrder = {}
        self.tOrder = []

    def empty(self):
        self.objectsByKey = {}
        self.O = {}
        self.oOrder = {}
        self.tOrder = []

    def __getitem__(self,key):
        if type(key) == tuple  and len(key) == 2:
            ot,oi = key
            return self.O[ot][oi]
        else:
            return self.oOrder[key]

    def getO(self, o_type,o_name):
        return self.O[o_type][o_name]

    def names(self, o_type):
        return [o.name for o in self.O[o_type].values()]

    def get_object_types(self):
        return self.tOrder
    
    def add(self,ot,oi,params=None,deps=None):
        o = OGO(ot,oi,params,deps)
        self.addObject(o)
        return o
 
    def addObject(self, ob):
        ob.OG = self
        key = str(ob.type)+":"+str(ob.name)

        if key in self.objectsByKey:
            raise ObjectGraphException(f"The object of type {ob.type} and key {ob.name} is already added\n");
        self.objectsByKey[key] = ob;
        if ob.type in self.O:
            self.oOrder[ob.type].append(ob)
            self.O[ob.type][ob.name] = ob
        else:
            self.O[ob.type] = {}
            self.tOrder.append(ob.type)
            self.O[ob.type][ob.name] = ob
            self.oOrder[ob.type] = [ob]
                     
    def save(self,outFile):
        f=open(outFile, "w") 

        f.write("{\n")            
        for ot in self.tOrder:
            for o in self.oOrder[ot]:
                f.write("\t\""+str(o.name)+"."+str(o.type)+ "\": {\n")
                f.write("\t\t\"id\": \""+str(o.name)+"\",\n")
                f.write("\t\t\"type\": \""+str(o.type)+"\",\n")
                f.write("\t\t\"deps\":["  +  ",".join(["\""+d.type+"/"+d.name+"\"" for d in o.deps]) + "],\n" )
                f.write("\t\t\"params\": {\n")
                out = ""
                for p,v in sorted(o.params.items()):
                    w = "\""+str(v)+"\""
                    out += "\t\t\t\""+p+"\": "+w+",\n"
                f.write(out[:-2]+"\n")
                f.write("\t\t}\n")
                f.write("\t},\n\n")
        f.write("\t\"dummy\": \"dummy\"\n")
        f.write("}\n")
        f.close()

    def print_stats(OG):
        print("Types:")
        for tp,tpOs in sorted(OG.O.items()):
            print("\t", tp, ":", len(tpOs)) 
    

def load_object_graph(fname):
    with open(fname) as f:
        OD = json.load(f, object_pairs_hook=OrderedDict)
        OG = ObjectGraph()
        for k in OD:
            if k == "dummy": break
            og = dict(OD[k])
            oi = og['id']
            ot = og['type']
            params = og['params'] if 'params' in og else []
            deps = og['deps'] if 'deps' in og else []
            dt = [d.split('/') for d in deps]
            OG.add(ot,oi, params, [OG[ot,oi] for ot,oi in dt])
        return OG

