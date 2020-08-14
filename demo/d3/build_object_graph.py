#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import sys

OG = ObjectGraph(".")

OG.addObject(OGO("base","o", {"a":"alabala nica"}))

for i in range(2):
    OG.addObject(OGO("sample","sample%d" % (i),{"g":"value%d" % (i)},[OG["base","o"]]))

OG.addObject(OGO("report","o", {}, OG["sample"]))

OG.execARGVcommands(sys.argv)

