add_targets('denovo_calls.vcf', 'igv-session.xml')

rule call_denovos:
  input:
    bams=DT("mdup.bam"),
    idx =DT("mdup.bam.bai"),
    ref = DT("chrAll.fa", dot="reference", level=3)
  output:
    T("denovo_calls.vcf")
  params:
    bed = PP("target")
  log:
    **(LFS("denovo_calls.vcf"))
  shell:
    "call_denovo.py {input.bams} {input.ref} {params.bed} {wildcards.oid} >{output} 2>{log.E}"

rule igv:
  input:
    bams=DT("mdup.bam"),
    idx =DT("mdup.bam.bai"),
    ref = DT("chrAll.fa", dot="reference", level=3) 
  output:
    T("igv-session.xml")
  log:
    **(LFS("igv-session.xml"))
  shell:
    " makeIGVsession.py trio `dirname {output}` {input.bams} 2>{log.E} && [ -s {output} ]"

rule clean_trio:
  shell:
    "cd trio; for n in *; do (cd $n; rm -rf *) done"
