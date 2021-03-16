add_targets("sample.bam",
        "sorted_reads.bam",
        "sorted_reads.bam.bai")

rule bwa_map:
    input:
        PP("ref"),
        PP('samplesDir') + "/{oid}.fastq"
    output:
        T("sample.bam")
    shell:
        "bwa mem {input} | samtools view -Sb - > {output}"

rule samtools_sort:
    input:
        T("sample.bam")
    output:
        T("sorted_reads.bam")
    shell:
        "samtools sort -T {input} -O bam {input} > {output}"

rule samtools_index:
    input:
        T("sorted_reads.bam")
    output:
        T("sorted_reads.bam.bai")
    shell:
        "samtools index {input}"
