rule l2_allB:
  input:
    lambda wc : DT(wc,'B.txt',"level1"),
    lambda wc : DT(wc,'a.txt',"base")
  log:  **(EFS('allB.txt'))
  output:
    T("allB.txt")
  shell:
    "(time echo {input} > {output} \
            2> {log.E} \
     ) 2> {log.T}"


rule level2_obj:
  input: 
        lambda wc: DT(wc,"obj.flag"),
        T('allB.txt')
  output:
    T("obj.flag")
  shell:
    "touch {output}"



