import pandas as pd
def run(proj, OG):
  PED = pd.read_table(proj.parameters['pedigree'], 
                      sep='\t', header=0)
  OG.add('reference','hg38', {"symlink.ref.fa":
                              proj.parameters["reference"], "symlink.ref.fa.fai":
                              proj.parameters["refIndex"]})
  for i, r in PED.iterrows():
    OG.add('individual',r['personId'], 
               {'fqId':r['fastqId']},OG['reference'])
  for i, r in PED.iterrows():
    if r['fatherId'] == '.' or r['motherId'] == '.': 
       continue
    OG.add('trio',r['personId'], {},
           [OG['individual',r[i]] 
            for i in ['fatherId', 'motherId','personId']])
  OG.add('denovo','all', {}, OG['trio'])
