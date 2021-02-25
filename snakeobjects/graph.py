#!/bin/env python

from snakeobjects import Project, ObjectGraph, load_object_graph
import random,sys

# project = Project("variants/batch3")
# OG = project.objectGraph
# OG = load_object_graph("objects/.snakeobjects/OG.json")
# OG = load_object_graph("OG-rnaSeq.json")


clrsStr = '''
red green blue yellow salmon orange purple mediumslateblue linen powderblue pink skyblue
wheat mintcream magenta transparent olivedrab yellowgreen thistle
navajowhite sienna tan purple saddlebrown 
palevioletred orchid steelblue slategray orangered mediumaquamarine seagreen
mediumpurple paleturquoise slateblue tomato limegreen palegoldenrod midnightblue
turquoise whitesmoke mediumturquoise white navyblue mediumorchid mediumvioletred
plum peru peachpuff oldlace seashell 
mediumseagreen mistyrose mediumblue sandybrown moccasin violet snow
rosybrown violetred royalblue maroon papayawhip mediumspringgreen springgreen
palegreen slategrey
'''

colors = clrsStr.split()
colors = [x for x in colors if x[-1] not in "012345689"]

def plotGraph(OG, width=0.10, penwidth=0.1, arrowsize=0.1):
    print("graph:",width, penwidth, arrowsize, file=sys.stderr)
    print("digraph digraphname {")
    print('''
    graph [ size = "60,60" ];

    node [shape = circle,
      style = filled,
    ''')
    print(
      "width = %f," % width
        )
    print('''
      color = grey,
      label = ""]
    ''')

    def o2K(o):
        return ('"' + o.oType + "_" + o.oId + '"').replace('.','_')

    assignedColors = set()
    for t,c in zip(OG.get_object_types(),colors):
    
        print(f"{t} -> {c}",file=sys.stderr)
        print(f"node [fillcolor = {c}]")
        for o in OG[t]:
            print(o2K(o))
        print()

    print('''
    edge [color = black,
    ''')
    print(
          "penwidth= %f," % penwidth,
          "arrowsize= %f," % arrowsize
    )
    print("]")

    for t in OG.tOrder:
        for o in OG[t]:
            for d in o.deps:
                print (o2K(d) + " -> " + o2K(o) + ";")

    print()
    print("subgraph cluster_1 {")
    print('''
           label=Legend
	   nodesep=0.02
	   ranksep="0.02 equally"

	   node [ label="", fontsize=10, width=.2, shape=circle] 
          ''')
    for t,c in zip(OG.get_object_types(),colors):
    
        print(f"node [ fillcolor=%s ] %s_s" % (c, t))
    print()

    print("node [shape=none, fillcolor=none]")
    print()

    for t in OG.get_object_types():
    
        print("node [ label=%s ] %s_n" % (t, t))

    print()
    for t in OG.get_object_types():
    
        print("{ rank=same; %s_n;%s_s }" % (t, t))

    print()
    print("edge [style=invis]")
    for t in OG.get_object_types():
    
        print(f"%s_n -> %s_s" % (t, t))
    print()
    names = OG.get_object_types()
    for i in range(len(names)-1):
    
        print(f"%s_s -> %s_s" % (names[i], names[i+1]))
    print()
    print("}")
    
    print("}")

