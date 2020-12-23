#!/usr/bin/env python

from snakeobjects import Project, ObjectGraph

proj = Project()

OG = ObjectGraph()

OG.add("B","o", {"a":"alabala nica"})

OG.add('P','1',{'name':'Peter','dob':"3/10/2000"},OG['B'])
OG.add('P','2',{'name':'Paul' ,'dob':"4/20/2001"},OG['B'])
OG.add('P','3',{'name':'Mary' ,'dob':"5/30/2002"},OG['B'])
OG.add('P','4',{'name':'John' ,'dob':"6/11/2003"},OG['B'])

OG.add('AP','o',{},OG['P'])

OG.add('F','1',{'state':'happy'},[OG['P','1'],OG['P','2']])
OG.add('F','2',{'state':'  sad'},[OG['P','3'],OG['P','4']])
 
OG.add('AF','o',{},OG['F'] + OG['B'])

proj.prepare(OG)

