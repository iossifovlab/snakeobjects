add_targets('denovo_calls.txt')

rule call_denovos:
  input:
    bams=DT("mdup.bam"),
    idx =DT("mdup.bam.bai"),
    ref = DT("chrAll.fa", dot="reference", level=3)
  output:
    T("denovo_calls.txt")
  params:
    bed = PP("target")
  log:
    **(LFS("denovo_calls.txt"))
  shell:
    "call_denovo.py {input.bams} {input.ref} {params.bed} {wildcards.oid} >{output} 2>{log.E}"


rule clean_trio:
  shell:
    "cd trio; for n in *; do (cd $n; rm -rf *) done"
