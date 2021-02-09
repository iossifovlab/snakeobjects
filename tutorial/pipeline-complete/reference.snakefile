rule reference:
  input:
    T("chrAll.fa.amb"),
    T("chrAll.fa.ann"),
    T("chrAll.fa.bwt"),
    T("chrAll.fa.fai"),
    T("chrAll.fa.pac"),
    T("chrAll.fa.sa")
  output:
    touch(T("obj.flag"))

rule bwa_index:
  input:
    T("chrAll.fa")
  output:
    T("chrAll.fa.amb"),
    T("chrAll.fa.ann"),
    T("chrAll.fa.bwt"),
    T("chrAll.fa.pac"),
    T("chrAll.fa.sa")
  log:
    **(LFS("chrAll.fa.amb"))
  shell:
    "bwa index {input} 2> {log.E}"

rule reference_index:
  input:
    T("chrAll.fa")
  output:
    T("chrAll.fa.fai")
  shell:
    "samtools faidx {input}"

