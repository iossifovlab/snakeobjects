add_targets("fastq.bam", "sample.bam", "sample.bam.bai")

rule align:
    input:
        refFile       = DT("chrAll.fa"),
        refFileBwaIdx = DT("chrAll.bwaIndex.flag"),
    output:
        T("fastq.bam")
    params:
        fqId   = P("fqId")
    resources: 
        mem_mb = 5*1024
    threads: 5
    log: **LFS('align')    
    shell: 
        "(time bwa mem -t {threads} -R '@RG\\tID:{wildcards.oid}\\tSM:{wildcards.oid}' {input.refFile} ../input/fastq/{params.fqId}_1.fqz ../input/fastq/{params.fqId}_2.fqz 2> {log.E} |  samtools view -Sb - > {output}) 2> {log.T}"


rule reorganizedBam:
    input: T("fastq.bam"), DT("chrAll.fa")
    output: T("sample.bam"), T("markDupStats.txt")
    log: **LFS('reorganize')
    shell: '''
        (time 
            samtools fixmate -m --reference {input[1]} -O bam {input[0]} - |
            samtools sort -T {input[0]} -O bam | 
            samtools markdup -T {input[0]} -O bam -s --reference {input[1]} - {output[0]} 2> {output[1]}
        ) 2> {log.T}
    '''

rule indexBam:
    input: T("sample.bam")
    output: T("sample.bam.bai")
    log: **LFS('merge')
    shell: "(time samtools index -b {input} {output} > {log.O} 2> {log.E} ) 2> {log.T}" 
