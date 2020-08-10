#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])

OG.addObject(OGO("base","o", {}))

OG.execARGVcommands(sys.argv)

