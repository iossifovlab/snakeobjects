add_targets("bwa.flag")

rule bwa:
  input: T("chrAll.fa")
  output:T("bwa.flag")
  conda: "env-bwa.yaml"
  log: **(LFS("bwa.txt"))
  shell: "bwa index {input} -a bwtsw > {log.O} 2> {log.E} && touch {output}"

