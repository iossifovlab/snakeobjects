add_targets('B1.txt', 'B2.txt')

rule l1_big:
  input: 
        DT("a{r}.txt")
  output:
    T('B{r}.txt')
  log:  **(LFS('B{r}.txt'))
  params: thea=DP('a')
  shell:
    "(time echo {params.thea} > {output}) 2> {log.T}"



