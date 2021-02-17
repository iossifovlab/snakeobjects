import pandas as pd
from pathlib import Path
from collections import defaultdict

def run(proj, OG):
    fastqDir = Path(proj.parameters['fastqDir'])

    # creating the reference object
    OG.add('reference','o',{"symlink.chrAll.fa": proj.parameters['ref']})

    # creating the fastq objects
    personFastqs = defaultdict(list)
    fastqs = pd.read_table(proj.parameters["fastqMap"], sep='\t', header=0)
    for i, r in fastqs.iterrows():
        fqId = ".".join([r['flowcell'],r['lane'],r['barcode']])
        pId = r['individual']
        fqO = OG.add('fastq', fqId, {
                    'R1':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R1.fastq.gz",
                    'R2':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R2.fastq.gz",
                    'sampleId': pId,
                    'rg':       f"@RG\\\\tID:{fqId}\\\\tSM:{pId}"
                 }, OG['reference'])
        personFastqs[pId].append(fqO)

    OG.add('fastqSummary','o',deps=OG['fastq'])

    # creating the sample objects
    ped = pd.read_table(proj.parameters["pedigree"], sep='\t', header=0)
    for i, r in ped.iterrows():
        OG.add('sample',r['personId'], {
                    'familyId': r['familyId'],
                    'gender':   r['sex'],
                    'affected': r['affected']
                }, personFastqs[r['personId']])

    OG.add('sampleSummary','o',deps=OG['sample'])

    # creating the trio objects
    for i, r in ped.iterrows():
        if r['fatherId'] == '.' or r['motherId'] == '.': continue
        OG.add('trio',r['familyId'], { }, [
                    OG['sample',r['fatherId']], 
                    OG['sample',r['motherId']], 
                    OG['sample',r['personId']] 
                ])

    OG.add('trioSummary','o',deps=OG['trio'])


