def run(proj,OG):
    with open(proj.parameters['who_sings_what']) as f:
        for l in f:
            name, song = l.strip().split('\t')
            OG.add('P',name,{'name':name,'song':song})
