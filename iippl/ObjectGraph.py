#!/bin/env python


import sys
import copy
import re
import os 

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

class ObjectGraphException(Exception):
    pass

class ObjectGraph:

    def __init__(self,baseDir=None):
        if not baseDir:
            try:
                baseDir = os.environ['PROJECT_DIR']
            except KeyError:
                baseDir = '.'
        self.req = {}
        self.objectsByKey = {}
        self.O = {}
        self.oOrder = {}
        self.tOrder = []
        self.baseDir = baseDir 
        self.params = {}
        try:
            P = open(self.baseDir + "/parameters.txt")
            for l in P:
                l = l.strip()
                if not l:
                    continue
                if l[0] == '#':
                    continue
                k,v = l.split('=',1)
                self.params[k] = v
            P.close()
        except:
            pass

    def __getitem__(self,key):
        if len(key) == 2:
            ot,oi = key
            return self.O[ot][oi]
        else:
            return self.oOrder[key]

    def getO(self, o_type,o_name):
        return self.O[o_type][o_name]

    def getAll(self, o_type):
        return self.O[o_type].values()

    def createGlblMakefiel(self): 
        GLBL = open(self.baseDir + '/glbl.makefile', 'w')
        GLBL.write('all: \\\n')
        for tp,tpOs in sorted(self.O.items()):
            for ogo in list(tpOs.values()):
                ogo.dir = str(self.baseDir) + "/objLinks/" + str(ogo.type) + "/" + str(ogo.name)
                GLBL.write('\t' + ogo.dir + '/obj.flag \\\n')
        
        GLBL.write('\n')
        for tp,tpOs in sorted(self.O.items()):
            for ogo in list(tpOs.values()):
                GLBL.write(ogo.dir + '/obj.flag:')
                for d in ogo.deps:
                   GLBL.write(' ' + d.dir + '/obj.flag')
                GLBL.write('\n')
                GLBL.write('\tcd `dirname $@`; qqq_job.sh\n\n')
        GLBL.close()

    def createDirs(self):
        for tp,tpOs in sorted(self.O.items()):
            for ogo in list(tpOs.values()):
                if not os.path.exists(ogo.dir + "/log"):
                    os.makedirs(ogo.dir + '/log')


                for k,v in list(ogo.params.items()):
                    if not k.startswith("symlink."):
                        continue
                    dst = ogo.dir + "/" + k[8:]
                    src = v
                    os.system("ln -sf %s %s" % (src,dst))
                    # os.symlink(src,dst)

    def createParamsFiles(self):
        for tp,tpOs in sorted(self.O.items()):
            for ogo in list(tpOs.values()):
                if not os.path.exists(ogo.dir + "/log"):
                    os.makedirs(ogo.dir + '/log')

                PF = open(ogo.dir + '/params.txt', 'w')
                for k in sorted(ogo.params):
                    v = ogo.params[k]
                    if k.startswith("symlink."):
                        continue
                    PF.write(k + '=' + str(v) + '\n')

                PF.write("deps=" + " ".join([dep.dir for dep in ogo.deps]) + "\n")
                PF.write("parents=" + " ".join([par.dir for par in ogo.parents]) + "\n")

                PF.close()
        
    def fillHash(self, line, h):
        p = re.compile('^(.*)=(.*)$')
        m = p.match(line)
        if m:
            h[m.group()[:m.group().index("=")]] = m.group()[m.group().index("=")+1:]

    def printH(self, h):
        for nn in list(h.keys()):
            print(nn+"="+h[nn])

    def get_object_types(self):
        return self.tOrder
    
    def add(self,ot,oi,params=None,deps=None):
        o = OGO(ot,oi,params,deps)
        self.addObject(o)
        return o
 
    def addObject(self, ob):
        ob.OG = self
        ob.dir = str(self.baseDir) + "/objLinks/" + str(ob.type) + "/" + str(ob.name)
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

    def objDir(self, o):
        patt = self.params["obj.dir.pattern"]
        NNN = o.NNN
        nodeObjDir=patt.replace("${NNN}", str(NNN))
        return nodeObjDir+"/"+o.type+"/"+o.name

    def writeObjectGraph(self,outFile):
        f=open(self.baseDir + "/" + outFile, "w") 
        f.write("OBJECT_GRAPH\n")
        params = self.params
        params["baseDir"] = self.baseDir
        for k in sorted(self.params):
            f.write(k+"="+params[k]+"\n")
        f.write("\n")
            
        # ORIGINAL WRONG: for o in sorted(self.objectsByKey.values()):
        # variatn 1 (GOOD): for o in sorted(self.objectsByKey.values(),key=lambda x:x.name):
        #for oName,o in sorted([(x.name,x) for x in self.objectsByKey.values()]):
        #for o in sorted(list(self.objectsByKey.values()),key=lambda x:x.name):
        for ot in self.tOrder:
            for o in self.oOrder[ot]:
                f.write("OBJECT\n")
                f.write("id="+str(o.name)+"\n")
                f.write("type="+str(o.type)+"\n")
                f.write("dir="+str(o.dir)+"\n")
                if len(o.parents) > 0 and isinstance(o.parents[0],(OGO)):
                    f.write("parents=" + " ".join([par.dir for par in o.parents]) + "\n")
                else:
                    f.write("parents=" + " ".join([par for par in o.parents]) + "\n")
                if len(o.deps) > 0 and isinstance(o.deps[0],(OGO)):
                    f.write("deps="  +  " ".join([dep.dir for dep in o.deps]) + "\n" )
                else:
                    f.write("deps="  +  " ".join([dep for dep in o.deps]) + "\n" )
                f.write("PARAMS\n")
                for p,v in sorted(o.params.items()):
                    f.write(p+"="+str(v)+"\n")
                f.write("\n")
            
        f.write("END\n")
        f.close()
                     
    def writeObjectGraphJson(self,outFile):
        f=open(self.baseDir + "/" + outFile, "w") 
        f.write("{\n\t\"parameters\": {\n")
        params = self.params
        params["baseDir"] = self.baseDir
        out = ""
        for k in sorted(self.params):
            v = params[k] if params[k].isdecimal() else "\""+params[k]+"\""
            out += "\t\t\""+k+"\":" + v +",\n"

        f.write(out[:-2]+"\n")
        f.write("\t},\n\n")
            
        # ORIGINAL WRONG: for o in sorted(self.objectsByKey.values()):
        # variatn 1 (GOOD): for o in sorted(self.objectsByKey.values(),key=lambda x:x.name):
        #for oName,o in sorted([(x.name,x) for x in self.objectsByKey.values()]):
        #for o in sorted(list(self.objectsByKey.values()),key=lambda x:x.name):
        #for o in sorted(list(self.objectsByKey.values()),key=lambda x:(x.type,x.name)):
        for ot in self.tOrder:
            for o in self.oOrder[ot]:
                f.write("\t\""+str(o.name)+"."+str(o.type)+ "\": {\n")
                f.write("\t\t\"id\": \""+str(o.name)+"\",\n")
                f.write("\t\t\"type\": \""+str(o.type)+"\",\n")
                f.write("\t\t\"dir\": \""+str(o.dir)+"\",\n")
                if len(o.parents) > 0 and isinstance(o.parents[0],(OGO)):
                    f.write("\t\t\"parents\":\"" + ",".join([par.dir for par in o.parents]) + "\",\n")
                else:
                    f.write("\t\t\"parents\":\"" + ",".join([par for par in o.parents]) + "\",\n")
                if len(o.deps) > 0 and isinstance(o.deps[0],(OGO)):
                    f.write("\t\t\"deps\":\""  +  ",".join([dep.dir for dep in o.deps]) + "\",\n" )
                    f.write("\t\t\"deps_local\":["  +  ",".join(["\""+dep.type+"/"+dep.name+"\"" for dep in o.deps]) + "],\n" )
                    f.write("\t\t\"deps_type\":\""  + "/".join([dep.dir for dep in o.deps][0].split("/")[:-1]) + "\",\n" )
                    f.write("\t\t\"deps_objs\":["  +  ",".join(["\""+dep.name+"\"" for dep in o.deps]) + "],\n" )                
                else:
                    f.write("\t\t\"deps\":\""  +  ",".join([dep for dep in o.deps]) + "\",\n" )
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

    def loadFile(self,fname):
        hf = open(fname)

        done = 0
        if hf.readline().rstrip() != 'OBJECT_GRAPH':
            print("Not an OBJECT_GRAPH")
            raise
        for line in  hf:
            line = line.rstrip()
            if line == "":
                break
            self.fillHash(line, self.params);

        self.baseDir = self.params["baseDir"];
        
        for line in hf:
            line = line.rstrip()
            if line == "":
                done = 0;
                parents = []
                if req["parents"] != "":
                        parents = req["parents"].split(" ")
                deps = []
                if req['deps'] != "":
                        deps = req["deps"].split(" ")
                name = req['id']
                type = req["type"]
                ob = OGO(type,name)
                ob.dir = req["dir"]
                ob.deps = deps
                ob.parents = parents
                ob.params = pars
                
                # ob = { "name":req["id"], "type":req["type"], "dir":req["dir"], \
                #        "parents":parents, "deps":deps, "pars":pars, "NNN":req["NNN"] }
                
                self.addObject(ob)
            elif line == "OBJECT":
                req = {}
                pars = {}
                ob = {}
                continue
            elif line == "PARAMS":
                done = 1
                continue
            elif line == "END":
                pass
                # print "END\n"
            elif not done :
                self.fillHash(line, req);
            else:
                self.fillHash(line, pars);
                #self.printH(self.pars)
        hf.close()

    def printStats(OG):
        print("baseDir:", OG.baseDir)
        print("Parameters:")
        for k,v in sorted(OG.params.items()):
            print("\t", k + "=" + v)
        print("Types:")
        for tp,tpOs in sorted(OG.O.items()):
            print("\t", tp, ":", len(tpOs)) 
    

    def execARGVcommands(OG, ARGV):
        if len(ARGV) < 2:
            OG.printStats()
            return

        cmd = ARGV[1]

        if cmd == "help":
            print("Usage: ObjectGraph.py stats | createDirs <baseDir> | saveAs <fn> |  saveAsJson <fn> | printObjectsByType <type> |loadFile <fn>")
        elif cmd == 'stats':
            OG.printStats()
        elif cmd == 'createDirs':
            OG.createDirs()
            OG.writeObjectGraphJson("OG.json")
        elif cmd == 'saveAs':
            if len(ARGV) < 3:
                print("missing output file name")
                exit()
            fn = ARGV[2]
            OG.writeObjectGraph(fn)
        elif cmd == 'saveAsJson':
            if len(ARGV) < 3:
                print("missing output file name")
                exit()
            fn = ARGV[2]
            OG.writeObjectGraphJson(fn)
        elif cmd == 'loadFile':
            if len(ARGV) < 3:
                print("missing OG file name")
                exit()
            fn = ARGV[2]
            OG.loadFile(fn)
        elif cmd == 'printObjectsByType':
            if len(ARGV) < 3:
                print("missing type")
                exit()
            print("\n".join(OG.O[ARGV[2]]))
        else:
            print("unkown command:", cmd)


if __name__ == "__main__":        
    OG = ObjectGraph(".")

    OG.addObject(OGO("mIndelBT","o"))

    OG.addObject(OGO("mSnvDAE","o"))

    OG.addObject(OGO("fDesc","o"))

    fmIds = ["F_" + x for x in "a b c d e f".split()]
    for fmId in fmIds:
        OG.addObject(OGO("mIndelNucFam",fmId,
                        {"symlink.mainNucFam.txt": "/alabalaba/" + fmId },
                        [OG.O["mIndelBT"]["o"]]))
        
    OG.addObject(OGO("mIndelDAE","o",deps=list(OG.O["mIndelNucFam"].values())))
    OG.addObject(OGO("mDAE","o",deps= [OG.O["mIndelDAE"]["o"], OG.O["mSnvDAE"]["o"]]))

    OG.execARGVcommands(sys.argv)


    # OG.createDirs()
    # home work for Boris
    # OG.saveFile("V1.OG")
    # OG2 = loadFile("V1.OG")
    # OG2.saveFile("V2.OG")
    # assert diff sais noting for diff V1.OG V2.OG

