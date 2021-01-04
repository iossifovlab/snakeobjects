def run(proj,OG):
    n1 = int(proj.parameters['n1'])

    OG.add("base","o")
    for i in range(n1):
        OG.add("level1",str(i),{},OG["base"])
    OG.add("level2","o", {}, OG["base"] + OG["level1"])
