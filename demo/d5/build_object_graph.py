#!/usr/bin/env python

"""
GP
T,DT  [filter type type, and levels and mode]
P,DP
ERS
"""

from iippl.ObjectGraph import ObjectGraph, OGO
import os,sys

OG = ObjectGraph(os.environ['PROJECT_DIR'])

OG.addObject(OGO("B","o", {"a":"alabala nica"}))

OG.addObject(OGO('P','1',{'name':'Peter','dob':"3/10/2000"},OG['B']))
OG.addObject(OGO('P','2',{'name':'Paul' ,'dob':"4/20/2001"},OG['B']))
OG.addObject(OGO('P','3',{'name':'Mary' ,'dob':"5/30/2002"},OG['B']))
OG.addObject(OGO('P','4',{'name':'John' ,'dob':"6/11/2003"},OG['B']))

OG.addObject(OGO('AP','o',{},OG['P']))

OG.addObject(OGO('F','1',{'state':'happy'},[OG['P','1'],OG['P','2']]))
OG.addObject(OGO('F','2',{'state':'  sad'},[OG['P','3'],OG['P','4']]))
 
OG.addObject(OGO('AF','o',{},OG['F'] + OG['B']))

OG.execARGVcommands(sys.argv)

