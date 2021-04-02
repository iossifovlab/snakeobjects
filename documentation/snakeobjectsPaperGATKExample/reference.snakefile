add_targets('ref.fa', 'index.flag', 'ref.fa.fai', 'ref.dict')
rule indexRef:
  input:  PP('reference')
  output: T('ref.fa'),touch(T('index.flag'))
  shell:  "ln -s {input} {output[0]} &&    \
           bwa index {output[0]} -a bwtsw"

rule faidx:
  input:  T('ref.fa')
  output: T('ref.fa.fai')
  shell:  "samtools faidx {input}"

rule dict:
  input:  T('ref.fa')
  output: T('ref.dict')
  shell:  "samtools dict {input} -o {output}"
