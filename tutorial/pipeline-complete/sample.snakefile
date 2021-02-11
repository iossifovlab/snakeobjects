rule sample:
  input:
    DT('obj.flag'),
    T("ns.bam"),
    T("fixmate.bam"),
    T("cs.bam"),
    T("mdup.bam"),
    T("mdup.bam.bai"),
    T("coverage.png")
  output:
    touch(T("obj.flag"))

rule merge:
  input:
    DT("sample_ns.bam")
  output:
    T("ns.bam")
  log:
    **(LFS("ns.bam"))
  shell:
    "samtools merge -n {output} {input} 2>{log.E}"

rule fixmate:
  input:
    T("ns.bam")
  output:
    T("fixmate.bam")
  params:
    ref=PP("ref")
  log:
    **(LFS("fixmate.bam"))
  shell:
    "samtools fixmate -m --reference {params.ref} -O bam {input} {output} 2>{log.E}"

rule sample_cs:
  input:
    T("fixmate.bam")
  output:
    T("cs.bam")
  log:
    **(LFS("cs.bam"))
  shell:
    "samtools sort -T {input} -O bam {input} > {output} 2>{log.E}"

rule markdup:
  input:
    T("cs.bam")
  output:
    T("mdup.bam")
  params:
    ref = PP("ref")
  log:
    **(LFS("mdup.bam"))
  shell:
    "samtools markdup -T {input} -O bam -s --reference {params.ref} {input} {output} 2>{log.E}"

rule sample_idx:
  input:
    T("mdup.bam")
  output:
    T("mdup.bam.bai")
  shell:
    "samtools index -b {input} {output}"

rule depth:
  input:
    bam=T("mdup.bam"),
    idx=T("mdup.bam.bai")    
  output:
    T("depth.txt")
  params:
    target=PP("target"),
    ref = PP("ref")
  shell:
    "samtools depth -b {params.target} -q 30 -Q 30 {input.bam} > {output}"

rule sample_plot:
  input:
    T("depth.txt")
  output:
    T("coverage.png")
  log:
    **(LFS("coverage.png"))
  shell:
    "$SO_PIPELINE/coveragePlot.py {wildcards.oid} {output} {input} 2>{log.E}"

rule clean_sample:
  shell:
    "cd sample; for n in *; do (cd $n; rm -rf *) done"
