#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph,OGO
import sys

OG = ObjectGraph(".")

OG.addObject(OGO("index","o",{}))

for i in range(5):
    OG.addObject(OGO("sample",str(i),{"bamFile":"/somewhere/sample%d.bam" % (i)},[OG.O['index']['o']]))

OG.addObject(OGO("report","o",{},list(OG.O['sample'].values()) + [OG.O['index']['o']]))

OG.execARGVcommands(sys.argv)
