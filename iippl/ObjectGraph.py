#!/usr/bin/env python


import sys
import copy
import re
import os 
import json
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
        if type(key) == tuple  and len(key) == 2:
            ot,oi = key
            return self.O[ot][oi]
        else:
            return self.oOrder[key]

    def getO(self, o_type,o_name):
        return self.O[o_type][o_name]

    def names(self, o_type):
        return [o.name for o in self.O[o_type].values()]

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
                    os.system("touch -h -r %s %s" % (src,dst))
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
        '''
        params = self.params
        params["baseDir"] = self.baseDir
        for k in sorted(self.params):
            f.write(k+"="+params[k]+"\n")
        '''
        f.write(self.baseDir + "\n")
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
                f.write("deps="  +  " ".join([f'{d.type}:{d.name}' for d in o.deps]) + "\n" )
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

    def createSnakefile(self,ot,sfile):
        with open(sfile, 'w') as f:
            f.write("rule "+ ot + ":\n")
            f.write("  input:\n")
            f.write("    DT('obj.flag')\n") 
            f.write("  output:\n") 
            f.write("    touch(T('obj.flag'))\n\n")
 
    def writeMainSnakefile(self):
        pdir=os.environ["PROJECT_DIR"]
        with open(pdir+'/objLinks/main.snakefile', 'w') as f:
            #local='localrules: '+','.join(['all_'+ o for o in self.tOrder])
            #f.write(local+"\n\n")
            f.write("include: \""+ pdir +"/header.snakefile\"\n\n")
            f.write("rule all_main:\n")
            f.write("  input:\n")
            f.write("    expand(\"{od}/obj.flag\", od=all_obj_dirs())\n\n")

            for ot in self.tOrder:
                sfile = pdir + "/" + ot + ".snakefile"
                if not os.path.exists(sfile):
                    self.createSnakefile(ot,sfile)
                f.write("include: \""+pdir+"/" + ot + ".snakefile\"\n\n")
                f.write("rule all_" + ot + ":\n")
                f.write("  input:\n")
                f.write("    expand(\"" + ot + 
                        "/{o}/obj.flag\", o=OG.names(\""+ ot + "\"))\n\n")

     
    def printStats(OG):
        print("baseDir:", OG.baseDir)
        print("Parameters:")
        for k,v in sorted(OG.params.items()):
            print("\t", k + "=" + v)
        print("Types:")
        for tp,tpOs in sorted(OG.O.items()):
            print("\t", tp, ":", len(tpOs)) 
    

    def execARGVcommands(OG, ARGV=None):
        if not ARGV:
            ARGV = sys.argv

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
            OG.writeMainSnakefile()
            os.system("mkdir -p log")
            OG.writeObjectGraphJson("OG.json")
            # OG.writeObjectGraph("OG.OG")
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

def load_object_graph_json(fname):
    with open(fname) as f:
        OD = json.load(f, object_pairs_hook=OrderedDict)
        baseDir = dict(OD['parameters'])['baseDir']
        OG = ObjectGraph(baseDir)
        OD.pop('parameters')
        for k in OD:
            if k == "dummy":
                break
            og = dict(OD[k])
            oi = og['id']
            ot = og['type']
            params = og['params'] if 'params' in og else []
            deps = og['deps'] if 'deps' in og else []
            dt = [d.split('/') for d in deps]
            OG.add(ot,oi, params, [OG[ot,oi] for ot,oi in dt])
        return OG

def load_object_graph_og(fname):
    def fillHash(line, h):
        p = re.compile('^(.*)=(.*)$')
        m = p.match(line)
        if m:
            h[m.group()[:m.group().index("=")]] = m.group()[m.group().index("=")+1:]
        else:
            raise

    with open(fname) as hf:
        if hf.readline().rstrip() != 'OBJECT_GRAPH':
            print("Not an OBJECT_GRAPH")
            raise
        baseDir = hf.readline().rstrip()
        if hf.readline().rstrip() != "": 
            raise

        OG = ObjectGraph(baseDir)
        done = 0
        for line in hf:
            line = line.rstrip()
            # print("LLLL",line)
            if line == "":
                deps = []
                if req["deps"] != "":
                    for dpk in req["deps"].split(" "):
                        t,n = dpk.split(":")
                        deps.append(OG[t,n])
                n = req['id']
                t = req["type"]
                OG.add(t,n,pars,deps)
                done = 0
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
            elif not done :
                fillHash(line, req);
            else:
                fillHash(line, pars);
    return OG

def load_object_graph(fname):
    print (fname,file=sys.stderr)
    fn = os.path.basename(fname)
    if fn.endswith(".OG"):
        return load_object_graph_og(fname)
    elif fn.endswith(".json"):
        return load_object_graph_json(fname)
    else:
        print("load_object_graph: wrong file", fname)
        exit()

if __name__ == "__main__":        

    def p(name,ol):
        print(name,",".join([f"{o.type}:{o.name}" for o in ol]))

    '''
    OG = ObjectGraph()
    a = OG.add("A","o")
    b = OG.add("B","1",deps=OG['A'])
    b = OG.add("B","2",deps=OG['A'])
    b = OG.add("B","3",deps=OG['A'])
    c = OG.add("C","o",deps=OG['B'])

    p("b",b.deepDeps())
    p("a",a.deepDeps())
    p("c1",c.deepDeps())

    p("c2e",c.deepDeps(level=2))
    p("c2l",c.deepDeps(level=2,mode='lessOrEqual'))
    p("c2lA",c.deepDeps(dot="A",level=2,mode='lessOrEqual'))
    p("c2lB",c.deepDeps(dot="B",level=2,mode='lessOrEqual'))

    p("c3e",c.deepDeps(level=3))
    p("c3l",c.deepDeps(level=3,mode='lessOrEqual'))
    p("c3lA",c.deepDeps(dot="A",level=3,mode='lessOrEqual'))
    p("c3lB",c.deepDeps(dot="B",level=3,mode='lessOrEqual'))
    '''


    '''
    OG = ObjectGraph()

    OG.add("B","o", {"a":"alabala nica"})

    OG.add('P','1',{'name':'Peter','dob':"3/10/2000"},OG['B'])
    OG.add('P','2',{'name':'Paul' ,'dob':"4/20/2001"},OG['B'])
    OG.add('P','3',{'name':'Mary' ,'dob':"5/30/2002"},OG['B'])
    OG.add('P','4',{'name':'John' ,'dob':"6/11/2003"},OG['B'])

    OG.add('AP','o',{},OG['P'])

    OG.add('F','1',{'state':'happy'},[OG['P','1'],OG['P','2']])
    OG.add('F','2',{'state':'  sad'},[OG['P','3'],OG['P','4']])

    OG.add('AF','o',{},OG['F'] + OG['B'])

    OG['P','1'].deepDeps()
    '''

    # OG.execARGVcommands()

    OG = ObjectGraph(".")
    OG.add("B","o", {"a":"alabala nica"})
    OG.add('P','1',{'name':'Peter','dob':"3/10/2000"},OG['B'])
    OG.add('P','2',{'name':'Paul' ,'dob':"4/20/2001"},OG['B'])
    OG.add('P','3',{'name':'Mary' ,'dob':"5/30/2002"},OG['B'])
    OG.add('P','4',{'name':'John' ,'dob':"6/11/2003"},OG['B'])

    OG.add('AP','o',{},OG['P'])

    OG.add('F','1',{'state':'happy'},[OG['P','1'],OG['P','2']])
    OG.add('F','2',{'state':'  sad'},[OG['P','3'],OG['P','4']])

    OG.add('AF','o',{},OG['F'] + OG['B'])

    OG.writeObjectGraph("V1.OG")
    OG2 = load_object_graph("V1.OG")
    OG2.writeObjectGraph("V2.OG")
    os.system('diff V1.OG V2.OG')
    # assert diff sais noting for diff V1.OG V2.OG

