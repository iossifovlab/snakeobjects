import pandas as pd
def run(project, OG):
  PED = pd.read_table(project.parameters['pedigree'], 
                      sep='\t', header=0)
  OG.add('reference','hg38')
  for i, r in PED.iterrows():
    OG.add('individual',r['personId'], 
               {'fqId':r['fastqId']},OG['reference'])
  for i, r in PED.iterrows():
    if r['fatherId'] == '.' or r['motherId'] == '.': 
       continue
    OG.add('trio',r['personId'], {},
           [OG['individual',r[i]] 
            for i in ['fatherId', 'motherId','personId']])
  OG.add('denovos','all', {}, OG['trio'])
