import pandas as pd
from pathlib import Path

def run(proj, OG):
    fastqDir = Path(proj.parameters['fastqDir'])
    fastqs = pd.read_table(proj.parameters["fastqsFile"], sep='\t', header=0)

    OG.add('reference','o', {'symlink.chrAll.fa':proj.parameters['chrAllFile']})

    for i, r in fastqs.iterrows():
        OG.add('fastq', 
                ".".join([r['flowcell'],r['lane'],r['barcode']]),
                {
                    'R1':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R1.fastq.gz",
                    'R2':       fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R2.fastq.gz",
                    'sampleId': r['individual']
                },
                OG['reference']
             )
    OG.add('fastqSummary','o',deps=OG['fastq'])
