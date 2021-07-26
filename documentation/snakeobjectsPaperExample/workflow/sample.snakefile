add_targets("sample.bam","sorted_reads.bam","sorted_reads.bam.bai")

rule bwa_map:
    input:
        PP("ref"),
        EF("[PP:samplesDir][P:fqId]")
    output:
        T("sample.bam")
    conda: "environment.yml"
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        T("sample.bam")
    output:
        T("sorted_reads.bam")
    conda: "environment.yml"	
    shell:
        "samtools sort -T {input} -O bam {input} > {output}"

rule samtools_index:
    input:
        T("sorted_reads.bam")
    output:
        T("sorted_reads.bam.bai")
    conda: "environment.yml"	
    shell:
        "samtools index {input}"
