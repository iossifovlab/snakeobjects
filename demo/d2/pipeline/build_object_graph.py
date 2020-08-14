#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])
n1 = int(OG.params['n1'])

OG.addObject(OGO("base","o", {}))

for i in range(n1):
    OG.addObject(OGO("level1",str(i),{},[OG["base","o"]]))

OG.addObject(OGO("level2","o", {}, [OG["base","o"]] + OG["level1"]))

OG.execARGVcommands(sys.argv)

