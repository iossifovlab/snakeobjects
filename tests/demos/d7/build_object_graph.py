#!/usr/bin/env python

from snakeobjects import Project, ObjectGraph

proj = Project()

OG = ObjectGraph()

OG.add("B","o", {"a":"alabala nica"})

OG.add('P','1')
OG.add('P','2')


proj.prepare(OG)
