#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph

OG = ObjectGraph()

OG.add("base","o", {"a":"alabala nica"})

for i in range(2):
    OG.add("level1",str(i),{"g":"value%d" % (i)},OG["base"])

OG.add("level2","o", {}, OG["base"] + OG["level1"] )

OG.add("level3", "o", {}, [OG["level2","o"]] )

OG.execARGVcommands()

