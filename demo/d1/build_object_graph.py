#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])

OG.addObject(OGO("base","o", {"a":"alabala nica"}))

for i in range(2):
    OG.addObject(OGO("level1",str(i),{"g":"value%d" % (i)},[OG["base","o"]]))

OG.addObject(OGO("level2","o", {}, [OG["base","o"]] + OG["level1"] ))

OG.execARGVcommands(sys.argv)
