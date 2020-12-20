#!/usr/bin/env python

from snakeobjects.ObjectGraph import ObjectGraph
import sys

n = int(sys.argv[1])
sys.argv.pop(1)

OG = ObjectGraph()

OG.add("B","o", {"a":"alabala nica"})

for i in range(n):
    OG.add('P',str(i))


OG.execARGVcommands()

