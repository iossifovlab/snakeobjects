add_targets("a1.txt", "a2.txt", "b.txt")
add_targets("t.txt" )

rule build_a:
  output:
     t = T("t.txt")
  log:  **(LFS('t.txt'))
  params: a=P('a')
  shell:
    "(time echo {params.a}  > {output.t} \
           2> {log.E} \
     ) 2> {log.T}"

rule build_ab:
  output:
    T("a1.txt"), T("a2.txt"), T("b.txt")
  log:  **(LFS('b.txt'))
  shell:
    "touch {output[0]}; "
    "touch {output[1]}; "
    "touch {output[2]}; "




