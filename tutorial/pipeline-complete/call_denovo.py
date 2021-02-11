#!/usr/bin/env python

import numpy
from walker import *
from collections import defaultdict

bi = {x:i for i,x in enumerate('ACGT')}
ib = {i:x for i,x in enumerate('ACGT')}

class GT:
        pass

def genotype(gr):
    print(vars(gr), file=sys.stderr)
    denovo = 0
    for cnt in gr.C:
        if cnt[0] == 0 and cnt[1] == 0 and cnt[2] >5:
                denovo += 1
    if denovo:
        print ('\t'.join([gr.famId, 
                          gr.chrom, 
                          str(gr.position),
                          gr.refA, 
                          gr.altA, 
                          '/'.join([' '.join([str(x) for x in y]) 
                                    for y in gr.C])]))

def walk(P, minMaxCnt, famId):
    n = 0
    chrom = None

    for puc in iter(P):
        if len(puc) == 0:
                break
        assert len(set([puc[i].reference_pos 
                        for i in range(len(puc))])) == 1, 'walker is not sincronized on positions'
        chrom = puc[0].reference_name
        position = puc[0].reference_pos

        global goodData

        goodData = [defaultdict(int) for i in range(3)]

        n += 1

        for i,p in enumerate(puc):

            for  pileupread in p.pileups:
                read = pileupread.alignment
                qpos = pileupread.query_position
                qname = pileupread.alignment.query_name
    
                if pileupread.is_del or pileupread.is_refskip:
                    continue

                goodData[i][read.query_sequence[qpos]] += 1

            if min([len(x) for x in goodData]) == 0:
                    continue

        #print >>sys.stderr, 'chrom, position', chrom, position
        refA = P.fa.fetch(chrom, start=position, end=position+1)
        if not refA in 'ACGT':
            continue       

        cnt = [[0 if y not in x else x[y] for y in 'ACGT'] 
               for x in goodData]
        #print >>sys.stderr, position, cnt
        coverage = [sum([x for x in y]) for y in cnt]

        trio_coverage = sum(coverage)
        ref_coverage = sum([x[bi[refA]] for x in cnt])
        if ref_coverage == trio_coverage:
            continue

        alt_counts = [sum([0 if i == bi[refA] else x[i] 
                      for x in cnt]) for i in range(4)]
        alt_max =  max(alt_counts)
        if alt_max <= 2:
            continue
        altA = ib[alt_counts.index(alt_max)]

        refId = bi[refA]
        altId = bi[altA]

        total = [sum(x) for x in cnt]
        R = [x[refId] for x in cnt]
        A = [x[altId] for x in cnt]
        O = [total[i] - R[i] - A[i] for i in range(len(R))]
        C = numpy.array([R,A,O])
        gr = GT()
        gr.chrom = chrom
        gr.position = position
        gr.C = C
        gr.refA = refA
        gr.altA = altA
        gr.A = A
        gr.R = R
        gr.O = O
        gr.famId = famId

        genotype(gr)

    P.close()

def main():
    if len(sys.argv) <7:
        print("Usage: variantCaller.py <dad.bam mom.bam child.bam ref.fa bed file famId>")
        exit(1)
              
    global P
    files = []
    files = sys.argv[1:4]
    refF  = sys.argv[4]
    bedFn = sys.argv[5]
    famId = sys.argv[6]

    print('\t'.join('familyId chrom position refBase altBase baseCnts'.split(' ')))
    with open(bedFn, 'r') as f:
        for l in f:
            ch, st, en = l.strip().split('\t')
            P = PU(files, ch, int(st), int(en), refF, 30, 20)
            walk(P, 5, famId)

if __name__ == "__main__":
	main()
