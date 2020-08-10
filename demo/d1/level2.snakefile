rule l2_allB:
  input:
    lambda wc : DT(wc,'B.txt',"level1"),
    lambda wc : DT(wc,'a.txt',"base")
  output:
    T("allB.txt")
  shell:
    "echo {input} > {output}"


rule level2_obj:
  input: 
        lambda wc: DT(wc,"obj.flag"),
        T('allB.txt')
  output:
    T("obj.flag")
  shell:
    "touch {output}"



