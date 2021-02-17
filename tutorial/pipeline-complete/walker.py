#!/usr/bin/env python

import numpy
import pysam, sys, os
import copy

minMapQ = 30
minGoodBaseQ = 20

P = None

class PU:
    def __init__(pu, files, chrom, start, stop, refF, minMapQ, minGoodBaseQ):
        pu.refF = None
        if  os.path.isfile(refF):
            pu.fa = pysam.FastaFile(filename=refF, filepath_index=refF+'.fai')
        else:
            print("reference file %s" % refF, "is not found", file=sys.stderr) 

        pu.files = files
        pu.a = [pysam.AlignmentFile(x, mode='rc',reference_filename=refF) for x in files]

        if chrom:
            #print("chrom, start, stop", chrom, start, stop, file=sys.stderr)
            pu.pu = [x.pileup(chrom, start, stop,
                              stepper='samtools',
                              truncate = True,
                              min_mapping_quality=minMapQ, 
                              ignore_orphans=False,
                              ignore_overlaps=True,
                              max_depth=100000,
                              min_base_quality=20)  
                     for x in pu.a] 

        else:
            #print("no contig, no region", file=sys.stderr)
            pu.pu = [x.pileup(stepper='samtools',
                              truncate = True,
                              min_mapping_quality=minMapQ, 
                              ignore_orphans=False,
                              ignore_overlaps=True,
                              max_depth=100000,
                              min_base_quality=20)
                     for x in pu.a]

    def __iter__(pu):
        return pu

    def __next__(pu):
        try:
            puc = [next(x) for x in pu.pu]
        except ValueError:
            return []
        chrom_pos = [(x.reference_id,x.pos) for x in puc]
        M = max(chrom_pos)
        m = min(chrom_pos)
        while m != M:
            idx = chrom_pos.index(M)
            id = [i for i,cp in enumerate(chrom_pos) if i != idx]
            for i in id:
                while True:
                    try:
                        puc[i] = next(pu.pu[i])
                    except ValueError:
                        return []
                    if (puc[i].reference_id,puc[i].pos) >= M:
                        break
            chrom_pos = [(x.reference_id,x.pos) for x in puc]
            M = max(chrom_pos)
            m = min(chrom_pos)

        return puc

    def close(pu):
        for x in pu.a:
            x.close()
