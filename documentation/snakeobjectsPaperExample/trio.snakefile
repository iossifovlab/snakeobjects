add_targets("denovo_calls.txt")

rule callDenovos:
  input:
    bams=DT("sample.bam"),
    idx=DT("sample.bam.bai")
  output: T("denovo_calls.txt")
  params:
    bed = PP("target")
  shell:
    "call_denovo.py {input.bams} {params.bed} > {output}"

