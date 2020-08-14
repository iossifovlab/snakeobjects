#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph,OGO
import sys

OG = ObjectGraph(".")

OG.addObject(OGO("index","o",{}))

for i in range(5):
    OG.addObject(OGO("sample",str(i),{"bamFile":"/somewhere/sample%d.bam" % (i)},[OG['index','o']]))

OG.addObject(OGO("report","o",{},OG['sample'] + [OG.getO('index','o')]))

OG.execARGVcommands(sys.argv)
