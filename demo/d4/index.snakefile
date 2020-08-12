rule build_ref:
  output:
    T("ref.fa")
  log: **(EFS('build_ref'))
  shell:
    "(time touch {output} > {log.O} 2> {log.E}) 2> {log.T}"

rule index_obj:
  input:
    DT('obj.flag'),
    T("ref.fa")
  output: 
    touch(T('obj.flag'))
