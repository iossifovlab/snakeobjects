#!/usr/bin/env python

import sys,os
import myvariant

mv = myvariant.MyVariantInfo()

if len(sys.argv) < 2:
    print ("Usage: annotate.py <input vcf file>")
    exit(1)

inf = sys.argv[1]

with open(inf) as f:
    data = f.readlines()

data = [l.strip('\n\r').split('\t') for l in data ]

for i,r in enumerate(data):
    if r[0].startswith('#'):
        print('\t'.join(r))
    else:
        break
data = data[i:]

vars = [x[0]+':g.'+str(x[1])+x[3]+'>'+x[4] for x in data]

eff = [mv.getvariant(x, fields='snpeff', assembly='hg38') for x in vars]

for x,y in zip(data,eff):
        if y == None:
            x[-1] = x[-1] + ';;;'
        else:
            x[-1] = ';'.join([x[-1]]+
                             [z+'='+y['snpeff']['ann'][0][z] 
                              if type(y['snpeff']['ann']) == list 
                              else z+'='+y['snpeff']['ann'][z]  
                          for z in 'effect gene_id genename'.split(' ')])
        print('\t'.join(x))


