import pandas as pd
from pathlib import Path

def run(proj, OG):
    fastqDir = Path(proj.parameters['fastqDir'])
    fastqs = pd.read_table(proj.parameters["fastqsFile"], sep='\t', header=0)

    OG.add('ref','o', {'symlink.chrAll.fa':proj.parameters['chrAllFile']})

    for i, r in fastqs.iterrows():
        fqId = ".".join([r['flowcell'],r['lane'],r['barcode']])
        OG.add('fastq', fqId, 
                {
                    'R1':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R1.fastq.gz",
                    'R2':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R2.fastq.gz",
                    'sampleId': r['individual'],
                    'rg':       f"@RG\\\\tID:{fqId}\\\\tSM:{r['individual']}"
                },
                OG['ref']
             )
    OG.add('fastqSummary','o',deps=OG['fastq'])
