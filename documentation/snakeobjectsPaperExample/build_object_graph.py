import pandas as pd
from pathlib import Path
from collections import defaultdict

def run(proj, OG):
    OG.add('reference','hg38', {'symlink.chrAll.fa':proj.parameters['chrAllFile']})

    ped = pd.read_table(proj.parameters["pedigree"], sep='\t', header=0)
    for i, r in ped.iterrows():
        OG.add('individual',r['personId'], {"fqId":r['fastqId']},OG['reference'])

    for i, r in ped.iterrows():
        if r['fatherId'] == '.' or r['motherId'] == '.': continue
        OG.add('trio',r['personId'], 
               deps=[OG['individual',r[i]] for i in ['fatherId', 'motherId','personId']])
    OG.add('trioSummary','o',deps=OG['trio'])
    
