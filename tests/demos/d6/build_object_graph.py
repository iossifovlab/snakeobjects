#!/usr/bin/env python

from snakeobjects import Project,ObjectGraph
import sys

proj = Project()

n = int(sys.argv[1])
sys.argv.pop(1)

OG = ObjectGraph()
OG.add("B","o", {"a":"alabala nica"})
for i in range(n):
    OG.add('P',str(i))


proj.prepare(OG)

