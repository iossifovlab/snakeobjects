add_targets("sample.bam", "sample.bam.bai","markDupStats.txt")

rule merge:
    input: DT("fastq.bam")
    output: temp(T("raw.bam"))
    conda: "env-bwa.yaml"
    log: **LFS('merge')
    shell: "(time samtools merge -n {output} {input} > {log.O} 2> {log.E} ) 2> {log.T}"

rule reorganized_bam:
    input: T("raw.bam"), DT("chrAll.fa",level=2)
    output: T("sample.bam"), T("markDupStats.txt")
    conda: "env-bwa.yaml"
    log: **LFS('reorganize')
    shell: '''
        (time 
            samtools fixmate -m --reference {input[1]} -O bam {input[0]} - |
            samtools sort -T {input[0]} -O bam | 
            samtools markdup -T {input[0]} -O bam -s --reference {input[1]} - {output[0]} 2> {output[1]}
        ) 2> {log.T}
    '''

rule sample_idx:
    input: T("sample.bam")
    output: T("sample.bam.bai")
    conda: "env-bwa.yaml"
    log: **LFS('merge')
    shell: "(time samtools index -b {input} {output} > {log.O} 2> {log.E} ) 2> {log.T}" 

