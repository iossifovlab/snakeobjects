add_targets("pairNumber.txt","fastq.bam")

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

rule align:
    input:
        refFile       = DT("chrAll.fa"),
        refFileBwaIdx = DT("chrAll.bwaIndex.flag"),
        R1File        = P("R1"),
        R2File        = P("R2")
    output:
        T("fastq.bam")
    params:
        rg = P('rg')
    conda: "env-bwa.yaml"
    resources: 
        mem_mb = 5*1024,
    threads: 5
    log: **LFS('align')    
    shell:
        "(time bwa mem -t {threads} -R '{params.rg}' {input.refFile} {input.R1File} {input.R2File} 2> {log.E} | samtools view -Sb - > {output}) 2> {log.T}"

