import pandas as pd
from pathlib import Path
from collections import defaultdict

def run(proj, OG):
    fastqDir = Path(proj.parameters['fastqDir'])
    fastqs = pd.read_table(proj.parameters["fastqsFile"], sep='\t', header=0)

    OG.add('reference','o', {'symlink.chrAll.fa':proj.parameters['chrAllFile']})

    for i, r in fastqs.iterrows():
        fqId = ".".join([r['flowcell'],r['lane'],r['barcode']])
        OG.add('fastq', 
                fqId,
                {
                    'R1':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R1.fastq.gz",
                    'R2':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R2.fastq.gz",
                    'sampleId': r['individual'],
                    'rg':       f"@RG\\\\tID:{fqId}\\\\tSM:{r['individual']}"
                },
                OG['reference']
             )
    OG.add('fastqSummary','o',deps=OG['fastq'])

    sampleFastqOs = defaultdict(list)
    for o in OG['fastq']:
        sampleFastqOs[o.params['sampleId']].append(o) 
    for smId,fqOs in sampleFastqOs.items():
        OG.add('sample',smId,deps=fqOs)

    OG.add('sampleSummary','o',deps=OG['sample'])

    ped = pd.read_table(proj.parameters["pedigree"], sep='\t', header=0)
    for i, r in ped.iterrows():
        if r['fatherId'] == '.' or r['motherId'] == '.': continue
        if not set([r['fatherId'],r['motherId'],r['personId']]).issubset(set(sampleFastqOs)): continue
        OG.add('trio',r['personId'], {"familyId":r['familyId'],"sex":r['sex'],"affected":r['affected'] }, [
                    OG['sample',r['fatherId']],
                    OG['sample',r['motherId']],
                    OG['sample',r['personId']]
                ])
    OG.add('trioSummary','o',deps=OG['trio'])
    
