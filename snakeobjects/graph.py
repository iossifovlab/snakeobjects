#!/bin/env python

from snakeobjects import Project, ObjectGraph, load_object_graph
import random,sys,argparse

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

def plotGraph(OG, width=0.75, penwidth=1, arrowsize=1, legend=0, out='graph', id=0, shape='circle'):
    O = open(out+'.dot', 'w')
    print("graph:",width, penwidth, arrowsize, legend, out, file=sys.stderr)
    print("digraph digraphname {", file=O)
    print('''
    graph [ size = "60,60" ];

    node [shape = circle,
      style = filled,
    ''', file=O)
    print(
      "width = %f," % width
        , file=O)
    print('''
      color = grey,
      label = ""]
    ''', file=O)

    def o2K(o):
        return ('"' + o.oType + "_" + o.oId + '"').replace('.','_')

    assignedColors = set()
    for t,c in zip(OG.get_object_types(),colors):
    
        print(f"{t} -> {c}",file=sys.stderr)
        print(f"node [fillcolor = {c}, shape = {shape}]", file=O)
        for o in OG[t]:
            if id == 0:
                print(o2K(o), file=O)
            elif id == 1:
                print(f"node [label=\"{o.oId}\"] {o2K(o)}", file=O)
            elif id == 2:
                print(f"node [label=\"{o.oType}/{o.oId}\"] {o2K(o)}", file=O)               
        print('', file=O)

    print('''
    edge [color = black,
    ''', file=O)
    print(
          f"penwidth = {penwidth},",
          f"arrowsize = {arrowsize},"
    , file=O)
    print("]", file=O)

    for t in OG.tOrder:
        for o in OG[t]:
            for d in o.deps:
                print (o2K(d) + " -> " + o2K(o) + ";", file=O)

    print('', file=O)
    print('}', file=O)

    O.close()
    
    if legend:
        O = open(out+'_legend.dot', 'w')
        print("digraph digraphname {", file=O)
        print('''
        graph [ size = "60,60" ];

        node [shape = circle,
        style = filled,
        ''', file=O)
        print(
            "width = {width),"
            , file=O)
        print('''
        color = grey,
        label = ""]
        ''', file=O)
        print('''
           labelloc="t"
           label=Legend
	   nodesep=0.02
	   ranksep="0.02 equally"

	   node [ label="", fontsize=10, width=.2, shape=circle] 
          ''', file=O)
        for t,c in zip(OG.get_object_types(),colors):
    
            print(f"node [ fillcolor=%s ] %s_s" % (c, t), file=O)
        print(' ', file=O)

        print("node [shape=none, fillcolor=none]", file=O)
        print(' ', file=O)

        for t in OG.get_object_types():
    
            print("node [ label=%s ] %s_n" % (t, t), file=O)

        print(' ', file=O)
        for t in OG.get_object_types():
    
            print("{ rank=same; %s_n;%s_s }" % (t, t), file=O)

        print('', file=O)
        print("edge [style=invis]", file=O)
        for t in OG.get_object_types():
    
            print(f"%s_n -> %s_s" % (t, t), file=O)
        print('', file=O)
        names = OG.get_object_types()
        for i in range(len(names)-1):
    
            print(f"%s_s -> %s_s" % (names[i], names[i+1]), file=O)
        print('', file=O)
        print("}", file=O)
        O.close()

def driver(OG, data):
    print(data, file=sys.stderr)
    parser = argparse.ArgumentParser(prog='graph.py')

    parser.add_argument("-w", "--width", dest="width", default=.75, type=float, action='store',
        metavar="width", help="width of node, default is 0.75" )

    parser.add_argument("-p", "--penwidth", dest="penwidth", default=1, type=float, action='store',
        metavar="penwidth", help="thickness of edges, default is 1.0" )

    parser.add_argument("-a", '--arrowsize', dest='arrowsize', default=1, type=float, metavar='arrowsize',
        help='multiplicative scale factor for arrowheads, default is 1.0' )

    parser.add_argument("-l", '--legend', dest='legend', default=False, type=bool, metavar='legend',
        help='if True generates legend in a separate file, default is False' )

    parser.add_argument("-o", '--out', dest='out', default='graph', type=str, metavar='out',
        help='prefix of output file, default is graph' )

    parser.add_argument("-i", '--id', dest='id', default=0, type=int, metavar='id',
        help='place id labels in nodes: 0 - no id, 1 - oid, 2 - oType/oId' )

    parser.add_argument("-s", '--shape', dest='shape', default='circle', type=str, metavar='shape',
        help='shape of the node, default is circle, for all shape names see https://www.graphviz.org/doc/info/shapes.html' )
    args = parser.parse_args(data[1:])

    width = args.width
    penwidth = args.penwidth
    arrowsize = args.arrowsize
    legend = args.legend
    out = args.out
    id = args.id
    shape = args.shape
    
    plotGraph(OG, width=width, penwidth=penwidth, arrowsize=arrowsize, legend=legend, out=out, id=id, shape=shape)




