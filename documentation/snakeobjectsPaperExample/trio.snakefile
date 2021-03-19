add_targets("denovo_calls.txt")
rule callDenovos:
  input: bams=DT("sample.bam"), bai=DT("sample.bam.bai")
  output: T("denovo_calls.txt")
  params: tR = PP("targetRegions")
  shell: "call_denovo.py {input.bams} {params.tR} > {output}"

