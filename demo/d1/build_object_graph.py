#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])

OG.addObject(OGO("base","o", {}))

for i in range(2):
    OG.addObject(OGO("level1",str(i),{},[OG.O["base"]["o"]]))

OG.execARGVcommands(sys.argv)

