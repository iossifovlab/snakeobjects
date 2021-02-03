#!/bin/env python

from snakeobjects import Project, ObjectGraph, load_object_graph
import random,sys

# project = Project("variants/batch3")
# OG = project.objectGraph
# OG = load_object_graph("objects/.snakeobjects/OG.json")
# OG = load_object_graph("OG-rnaSeq.json")


clrsStr = '''
lightyellow3    lightyellow4    limegreen   linen   magenta
magenta1    magenta2    magenta3    magenta4    maroon
maroon1 maroon2 maroon3 maroon4 mediumaquamarine
mediumblue  mediumorchid    mediumorchid1   mediumorchid2   mediumorchid3
mediumorchid4   mediumpurple    mediumpurple1   mediumpurple2   mediumpurple3
mediumpurple4   mediumseagreen  mediumslateblue mediumspringgreen   mediumturquoise
mediumvioletred midnightblue    mintcream   mistyrose   mistyrose1
mistyrose2  mistyrose3  mistyrose4  moccasin    navajowhite
navajowhite1    navajowhite2    navajowhite3    navajowhite4       navy   
navyblue       none     oldlace olivedrab   olivedrab1
olivedrab2  olivedrab3  olivedrab4  orange  orange1
orange2 orange3 orange4 orangered   orangered1
orangered2  orangered3  orangered4  orchid  orchid1
orchid2 orchid3 orchid4 palegoldenrod   palegreen
palegreen1  palegreen2  palegreen3  palegreen4  paleturquoise
paleturquoise1  paleturquoise2  paleturquoise3  paleturquoise4  palevioletred
palevioletred1  palevioletred2  palevioletred3  palevioletred4  papayawhip
peachpuff   peachpuff1  peachpuff2  peachpuff3  peachpuff4
   peru        pink     pink1   pink2   pink3
pink4      plum     plum1   plum2   plum3
plum4   powderblue  purple  purple1 purple2
purple3 purple4    red         red1        red2   
   red3        red4     rosybrown   rosybrown1  rosybrown2
rosybrown3  rosybrown4  royalblue   royalblue1  royalblue2
royalblue3  royalblue4  saddlebrown salmon  salmon1
salmon2 salmon3 salmon4 sandybrown  seagreen
seagreen1   seagreen2   seagreen3   seagreen4   seashell
seashell1   seashell2   seashell3   seashell4   sienna
sienna1 sienna2 sienna3 sienna4 skyblue
skyblue1    skyblue2    skyblue3    skyblue4    slateblue
slateblue1  slateblue2  slateblue3  slateblue4  slategray
slategray1  slategray2  slategray3  slategray4  slategrey
   snow     snow1   snow2   snow3   snow4
springgreen springgreen1    springgreen2    springgreen3    springgreen4
steelblue   steelblue1  steelblue2  steelblue3  steelblue4
   tan         tan1        tan2        tan3        tan4   
thistle thistle1    thistle2    thistle3    thistle4
tomato  tomato1 tomato2 tomato3 tomato4
transparent turquoise   turquoise1  turquoise2  turquoise3
turquoise4  violet  violetred   violetred1  violetred2
violetred3  violetred4  wheat   wheat1  wheat2
wheat3  wheat4  white   whitesmoke  yellow
yellow1 yellow2 yellow3 yellow4 yellowgreen
'''

colors = clrsStr.split()
colors = [x for x in colors if x[-1] not in "012345689"]

def plotGraph(OG, width=0.05, penwidth=0.1, arrowsize=0.1):
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
    for t in OG.get_object_types():
        while True:
            c = random.choice(colors)
            if c in assignedColors: continue 
            assignedColors.add(c)
            break
    
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
    print("}")

