add_targets("quals.svg","all.vcf")

rule bcftools_call:
    input:
        fa=PP("ref"),
        bam=DT("sorted_reads.bam"),
        bai=DT("sorted_reads.bam.bai")
    output:
        T("all.vcf")
    shell:
        "samtools mpileup -g -f {input.fa} {input.bam} | bcftools call -mv - > {output}"

rule plot_quals:
    input:
        T("all.vcf")
    output:
        T("quals.svg")
    script:
        "plot-quals.py"
