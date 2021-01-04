def run(proj,OG,*args):
    n = int(args[0])

    OG.add("B","o", {"a":"alabala nica"})
    for i in range(n):
        OG.add('P',str(i))

