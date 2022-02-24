add_targets('index.flag')
rule indexRef:
  input:  T("ref.fa")
  output: touch(T('index.flag'))
  shell:  "bwa index {input} -a bwtsw"
