rule fastq:
  input:
    T("readNumber.txt"),
    T("fastq.bam"),
    DT("obj.flag")
  output:
    touch(T("obj.flag"))

rule align:
    input:
        refFile     = DT("chrAll.fa"),
        refFileBwaI = DT("chrAll.bwaIndex.flag"),
        R1File      = P("R1"),
        R2File      = P("R2")
    output:
        T("fastq.bam")
    params:
        rg = P('rg')
    conda: "env-bwa.yaml"
    shell:
        "bwa mem -R '{params.rg}' {input.refFile} {input.R1File} {input.R2File} | samtools view -Sb - > {output}"
 

rule countReads:
    input:
        P('R1'), P('R2')
    output:
        T("readNumber.txt")
    run:
        import gzip

        nReads = 0
        readLengths = set() 
        buff = []
        with gzip.open(input[0]) as R1F, \
             gzip.open(input[1]) as R2F:
            for l1,l2 in zip(R1F,R2F):
                buff.append((l1,l2))
                if len(buff) == 4:

                    # the read ids should be the same for read1 and read2
                    if buff[0][0].split()[0] != buff[0][1].split()[0]:
                        nReads = -2
                        break

                    # the sequences and the qulities for the two reads whould have the same length 
                    if len(buff[1][0]) != len(buff[3][0]) or \
                       len(buff[1][0]) != len(buff[1][1]) or \
                       len(buff[1][0]) != len(buff[3][0]):
                        nReads = -3
                        break
        
                    readLengths.add(len(buff[1][0]))            
                    nReads += 1
                    buff = []
            # if the number of lines in the fastq files is not devisible by 4
            if buff:
                nReads = -4

            # if there are reads of different length in the fastq files.
            if len(readLengths) != 1:
                nReads = -5

            readLength = ",".join(map(str,readLengths))
        
        with open(output[0],"w") as OF:
            OF.write(f'{wildcards.oid}\t{nReads}\t{readLength}\n')
