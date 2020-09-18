rule l1_big:
  input: 
        DT("a{r}.txt")
  output:
    T('B{r}.txt')
  log:  **(EFS('B{r}.txt'))
  params: thea=DP('a')
  shell:
    "(time echo {params.thea} > {output}) 2> {log.T}"

rule level1_obj:
  input: 
    DT("obj.flag"),
    T('B1.txt'),
    T('B2.txt')
  output:
    touch(T("obj.flag"))



