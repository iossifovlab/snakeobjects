#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph
import sys

OG = ObjectGraph()

OG.add("B","o", {"a":"alabala nica"})

OG.add('P','1')
OG.add('P','2')


OG.execARGVcommands()

