add_targets("sample.bam", "sample.bam.bai","markDupStats.txt")

rule merge:
    input: DT("fastq.bam")
    output: T("raw.bam")
    conda: "env-bwa.yaml"
    shell: "samtools merge -n {output} {input}"

rule reorganized_bam:
    input: T("raw.bam"), DT("chrAll.fa",level=2)
    output: T("sample.bam"), T("markDupStats.txt")
    conda: "env-bwa.yaml"
    shell: 
        '''
        samtools fixmate -m --reference {input[1]} -O bam {input[0]} - |
        samtools sort -T {input[0]} -O bam | 
        samtools markdup -T {input[0]} -O bam -s --reference {input[1]} - {output[0]} 2> {output[1]}
        '''

rule sample_idx:
    input: T("sample.bam")
    output: T("sample.bam.bai")
    conda: "env-bwa.yaml"
    shell: "samtools index -b {input} {output}"

