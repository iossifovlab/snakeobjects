rule l2_allB:
  input:
    DT('B.txt',"level1"),
    DT('a.txt',"base")
  log:  **(EFS('allB.txt'))
  output:
    T("allB.txt")
  shell:
    "(time echo {input} > {output} \
            2> {log.E} \
     ) 2> {log.T}"

rule level2_obj:
  input: 
    DT("obj.flag"),
    T('allB.txt')
  output:
    touch(T("obj.flag"))

