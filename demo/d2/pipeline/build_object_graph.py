#!/usr/bin/env python

from iippl.ObjectGraph import ObjectGraph
import yaml,os

CF = open(os.envirno["PROJECT_DIR"] + "/parameters.yaml", 'r')
config = yaml.safe_load(CF)
CF.close()

OG = ObjectGraph()

n1 = int(config['n1'])

OG.add("base","o")

for i in range(n1):
    OG.add("level1",str(i),{},OG["base"])

OG.add("level2","o", {}, OG["base"] + OG["level1"])

OG.execARGVcommands()

