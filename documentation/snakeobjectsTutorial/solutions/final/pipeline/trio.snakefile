add_targets("denovo_calls.txt")

rule callDenovos:
  input:
    bams=DT("sample.bam"),
    idx=DT("sample.bam.bai")
  output: T("denovo_calls.txt")
  conda: "env-pysam.yaml"
  params:
    bed = PP("target")
  log: **LFS("denovo_calls")
  shell:
    "(time call_denovo.py {input.bams} {params.bed} > {output} 2>{log.E}) 2> {log.T}"

