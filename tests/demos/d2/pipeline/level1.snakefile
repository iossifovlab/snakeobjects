add_targets('B.txt')

rule l1_big:
  input: 
        DT("a.txt")
  output:
    T('B.txt')
  log:  **(LFS('B.txt'))
  shell:
    "(time touch {output}) 2> {log.T}"


