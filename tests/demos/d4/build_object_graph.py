#!/usr/bin/env python

from snakeobjects import Project, ObjectGraph

proj = Project()

OG = ObjectGraph()
OG.add("index","o",{})
for i in range(5):
    OG.add("sample",str(i),{"bamFile":"/somewhere/sample%d.bam" % (i)},OG['index'])
OG.add("report","o",{},OG['sample'] + OG['index'])

proj.prepare(OG)
