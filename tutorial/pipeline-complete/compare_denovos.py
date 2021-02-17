#!/usr/bin/env python
from collections import defaultdict

files = ['../input/denovos.txt', '../project-complete/objects/trioSummary/o/denovo_calls.vcf']

denovos = [defaultdict(),defaultdict()]

for i,fn in enumerate(files):
    with open (fn) as f:
        for l in f:
            if l.startswith('#'): continue
            cs = l.strip().split('\t')
            if '->' in l:
                refA = l[l.index('->')-1]
                altA = l[l.index('->')+2]
                
                denovos[i][cs[1]] = [refA,altA, cs[3]]
            else:
                denovos[i][cs[1]] = cs[2:]

common = set.intersection(set(denovos[0].keys()),set(denovos[1].keys()))
for k in common:
    assert( denovos[0][k][:2] == denovos[1][k][1:3]) 

not_found = set.difference(set(denovos[0].keys()),set(denovos[1].keys()))
false_pos = set.difference(set(denovos[1].keys()),set(denovos[0].keys()))

print ('not_found:')
for k in not_found:
    print('\t'.join([k]+denovos[0][k]))
print ('false positives:')
for k in false_pos:
    print('\t'.join([k]+denovos[1][k]))


