add_targets("chrAll.fa.amb", \
            "chrAll.fa.ann", \
            "chrAll.fa.bwt", \
            "chrAll.fa.fai", \
            "chrAll.fa.pac", \
            "chrAll.fa.sa", \
            "chrAll.fa.fai")

rule bwa_index:
  input:
    T("chrAll.fa")
  output:
    T("chrAll.fa.amb"),
    T("chrAll.fa.ann"),
    T("chrAll.fa.bwt"),
    T("chrAll.fa.pac"),
    T("chrAll.fa.sa")
  conda:
    "envs/bwa.yaml"
  log:
    **(LFS("chrAll.fa.amb"))
  shell:
    "bwa index {input} 2> {log.E}"

rule reference_index:
  input:
    T("chrAll.fa")
  output:
    T("chrAll.fa.fai")
  conda:
    "envs/bwa.yaml"
  shell:
    "samtools faidx {input}"

