
rule build_a:
  output:
     t = T("t.txt")
  log:  **(EFS('t.txt'))
  shell:
    "(time touch {output.t} \
           > {log.O}  \
           2> {log.E} \
     ) 2> {log.T}"

rule build_ab:
  output:
    T("a.txt"), T("b.txt")
  log:  **(EFS('b.txt'))
  shell:
    "touch {output[0]}; "
    "touch {output[1]}; "

rule base_obj:
  input:
    T("a.txt"), T("b.txt"), T("t.txt")
  output:
    touch(T("obj.flag"))



