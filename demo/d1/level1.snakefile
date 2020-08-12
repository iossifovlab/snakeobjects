rule l1_big:
  input: 
        DT("a.txt")
  output:
    T('B.txt')
  log:  **(EFS('B.txt'))
  params: thea=DP('a')
  shell:
    "(time echo {params.thea} > {output}) 2> {log.T}"

rule level1_obj:
  input: 
    DT("obj.flag"),
    T('B.txt')
  output:
    touch(T("obj.flag"))



