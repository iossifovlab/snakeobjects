rule l1_big:
  input: 
        lambda wc: DT(wc,"a.txt")
  output:
    T('B.txt')
  shell:
    "touch {output}"

rule level1_obj:
  input: 
        lambda wc: DT(wc,"obj.flag"),
        T('B.txt')
  output:
    T("obj.flag")
  shell:
    "touch {output}"



