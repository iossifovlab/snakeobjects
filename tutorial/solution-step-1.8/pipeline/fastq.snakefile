add_targets("pairNumber.txt")

rule countReads:
    input: P('R1'), P('R2')
    output: T("pairNumber.txt")
    run:
        import gzip

        nPairs= 0
        buff = []
        with gzip.open(input[0]) as R1F, \
             gzip.open(input[1]) as R2F:
            for l1,l2 in zip(R1F,R2F):
                buff.append((l1,l2))
                if len(buff) == 4:
                    nPairs += 1
                    buff = []
        
        with open(output[0],"w") as OF:
            OF.write(f'{wildcards.oid}\t{nPairs}\n')
