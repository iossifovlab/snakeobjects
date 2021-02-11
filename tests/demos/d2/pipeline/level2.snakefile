add_targets('allB.txt')

rule l2_allB:
  input:
    DT('B.txt',"level1"),
    DT('a.txt',"base")
  log:  **(LFS('allB.txt'))
  output:
    T("allB.txt")
  shell:
    "(time echo {input} > {output} \
            2> {log.E} \
     ) 2> {log.T}"




