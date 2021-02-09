rule trio:
  input:
    DT("obj.flag"),
    T("denovo_calls.txt")
  output:
    touch(T("obj.flag"))


rule call_denovos:
  input:
    DT("mdup.bam")
  output:
    T("denovo_calls.txt")
  params:
    ref = PP("ref"),
    bed = PP("target")
  log:
    **(LFS("denovo_calls.txt"))
  shell:
    "$SO_PIPELINE/call_denovo.py {input} {params.ref} {params.bed} {wildcards.oid} >{output} 2>{log.E}"


rule clean_trio:
  shell:
    "cd trio; for n in *; do (cd $n; rm -rf *) done"
