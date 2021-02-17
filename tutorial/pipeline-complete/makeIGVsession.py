#!/bin/env python
import sys
import os
import numpy as np
import subprocess
from string import Template
from glob import glob
import yaml

ref = 'hg38'

if len(sys.argv) <6:
    print("Usage: makeIGVsession.py <object type> <directory> <bam files>", file=sys.stderr)
    exit(1)
else:
    objType = sys.argv[1]
    objDir = sys.argv[2]
    dad, mom, child = sys.argv[3:]

if os.path.exists(objDir):
    os.chdir(objDir)
else:
    print("path does not exist:", objDir, file=sys.stderr)
    exit(1)

start_session='<?xml version="1.0" encoding="UTF-8" standalone="no"?> \
<Session genome="' + ref + '" hasGeneTrack="true" hasSequenceTrack="true" locus="chr1:144814725-144814783" version="7">'

resource=Template('        <Resource path="$NAME"/>')

start_panel=Template('    <Panel height="$HIGHT" name="Panel$NUMBER" width="$WIDTH">')

track1=Template('        <Track altColor="0,0,178" autoScale="true" color="175,175,175" colorScale="ContinuousColorScale;0.0;302.0;255,255,255;175,175,175" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="$PATH_COVERAGE" name="$NAME Coverage" showReference="false" snpThreshold="0.2" sortable="true" visible="true"> \
            <DataRange baseline="0.0" drawBaseline="false" flipAxis="false" maximum="302.0" minimum="0.0" type="LINEAR"/>\n\
        </Track>')

track2=Template('        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="EXPANDED" featureVisibilityWindow="-1" fontSize="10" id="$PATH" name="$NAME" showSpliceJunctions="false" sortable="true" visible="true"> \
            <RenderOptions colorByTag="" colorOption="SAMPLE" flagUnmappedPairs="false" groupByTag="" maxInsertSize="1000" minInsertSize="50" shadeBasesOption="QUALITY" shadeCenters="true" showAllBases="false" sortByTag=""/> \
        </Track>\n\
    </Panel>')

feature_panel='    <Panel height="118" name="FeaturePanel" width="1574"> \
        <Track altColor="0,0,178" autoScale="false" color="0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" id="Reference sequence" name="Reference sequence" sortable="false" visible="true"/> \
        <Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;308.0;255,255,255;0,0,178" displayMode="COLLAPSED" featureVisibilityWindow="-1" fontSize="10" height="35" id="' + ref + '_genes" name="RefSeq Genes" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"> \
            <DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="308.0" minimum="0.0" type="LINEAR"/> \
        </Track> \
    </Panel>'

layout_panel=Template('    <PanelLayout dividerFractions="$DIVIDERS"/> \
    <HiddenAttributes> \
        <Attribute name="NAME"/> \
        <Attribute name="DATA FILE"/> \
        <Attribute name="DATA TYPE"/> \
    </HiddenAttributes>')

#path='C:\\Users\\boris\Desktop\New folder\ '

def makeResources(names, path, f):
    global start_session
    f.write(start_session)
    f.write('    <Resources>\n')
    for n,p in zip(names, path):
        tmp = p
        f.write(resource.substitute(NAME=tmp) + '\n')
    f.write('    </Resources>\n')
    return

def makePanel(i, n, path, f):
    global start_panel
    f.write(start_panel.substitute(HIGHT='3000', NUMBER=str(i+1), WIDTH='1500') + '\n')
    tmp = path[i] + '_coverage'
    f.write(track1.substitute(PATH_COVERAGE=tmp, NAME=n) + '\n')
    tmp = path[i]
    f.write(track2.substitute(PATH=path[i], NAME=path[i]) + '\n')
    return

def makeSession(dirn):
    global feature_panel
    global layout_panel
    global objDir
    global objType


    names = ['father', 'mother', 'child']
    path = ['../../'+dad, '../../'+mom, '../../'+child]
    #path = subprocess.check_output('pwd', shell=True, stderr=subprocess.STDOUT).strip()
    f = open('igv-session.xml', 'w')
    makeResources(names, path, f)
    for i,n in enumerate(names):
        makePanel(i,n,path, f)
    f.write(feature_panel + '\n')
    a = .90/(len(names) + 1)
    d = ','.join(map(str, [0.005 + x*a for x in range(len(names) +1)]))
    f.write(layout_panel.substitute(DIVIDERS=d) + '\n')
    f.write('</Session>\n')
    f.close()

    return names

makeSession(objDir)
    
'''
i=1
for obj in os.listdir('.'):
    print i
    i = i+1
    makeSession(obj)
'''

