#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])
n1 = int(OG.params['n1'])

OG.addObject(OGO("base","o", {}))

for i in range(n1):
    OG.addObject(OGO("level1",str(i),{},[OG.O["base"]["o"]]))

OG.addObject(OGO("level2","o", {}, [OG.O["base"]["o"]] + list(OG.O["level1"].values())))

OG.execARGVcommands(sys.argv)

