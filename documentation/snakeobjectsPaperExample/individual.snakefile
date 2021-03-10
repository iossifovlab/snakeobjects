add_targets("fastq.bam", "sample.bam", "sample.bam.bai")

rule align:
    input: ref=DT("chrAll.fa"), refIdx=DT("chrAll.bwaIndex.flag")
    output: T("fastq.bam")
    params:
            fD=PP('fastqDir'), fI=P('fqId')
    shell: 
        "bwa mem -R '@RG\\tID:{wildcards.oid}\\tSM:{wildcards.oid}' {input.ref} {params.fD}/{params.fI}_1.fqz {params.fD}/{params.fI}_2.fqz | \
         samtools view -Sb - > {output}"

rule reorganizedBam:
    input: T("fastq.bam")
    output: T("sample.bam")
    shell: "samtools sort {input} -O bam > {output}"

rule indexBam:
    input: T("sample.bam")
    output: T("sample.bam.bai")
    shell: "samtools index -b {input} {output}"
