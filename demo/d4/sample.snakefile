rule realign:
  input:
    DT('ref.fa')
  output:
    T('realigned.bam')
  params:
      bam=P('bamFile')
  log: **EFS('realing')
  shell:
    "(time echo realinged {params.bam} using {input} > {output} 2> {log.E}) 2> {log.T}"

rule sample_obj:
  input:
    DT('obj.flag'),
    T('realigned.bam')
  output: 
    touch(T('obj.flag'))
