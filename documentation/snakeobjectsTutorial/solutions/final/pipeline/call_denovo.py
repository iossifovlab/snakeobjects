#!/usr/bin/env python

import pysam,sys
import pandas as pd
import numpy as np
from collections import Counter,defaultdict

dadBamF,momBamF,chlBamF,targetFile = sys.argv[1:]

AFS = [pysam.AlignmentFile(bmf) for bmf in [dadBamF,momBamF,chlBamF]]
childId ,= {x['SM'] for x in AFS[2].header['RG']}

print("\t".join(['childId','chrom','pos','newAllele'] + [f'{p}.cnt.{a}' for p in ['mom','dad','chl'] for a in "ACGT"]))
TT = pd.read_table(targetFile, sep='\t',names=['chr','beg','end'])
for ri,rgn in TT.iterrows():
    cntBuff = defaultdict(list)
    for AFi,AF in enumerate(AFS):
        plps = AF.pileup(rgn['chr'],rgn['beg'],rgn['end'])
        for plp in plps:
            cnt = Counter([n.upper() for n in plp.get_query_sequences()])
            cntA = np.array([cnt[n]  for n in 'ACGT'])
            cntBuff[plp.reference_pos].append(cntA)
    for pos,cntAs in cntBuff.items():
        if len(cntAs) != 3: continue
        trioCnt = np.vstack(cntAs)
        dpths = trioCnt.sum(axis=1)
        for dnvAlleleCandidate in range(4):
            if trioCnt[2,dnvAlleleCandidate] < 3: continue # the allele should be seen 3 or more time in the child
            if trioCnt[0,dnvAlleleCandidate] > 0: continue # the allele should NOT be seen in the parents
            if trioCnt[1,dnvAlleleCandidate] > 0: continue # the allele should NOT be seen in the parents
            if dpths[0] < 10 or dpths[1] < 10: continue    # each of the parents should have at least 10 reads at the position 
            print("\t".join(map(str,[childId, rgn['chr'], pos, "ACGT"[dnvAlleleCandidate]] + list(trioCnt.flatten()))))
for AF in AFS: AF.close()
