rule l2_allB:
  input:
    DT('B.txt',"level1"),
    DT('a.txt',"base")
  params: 
        allGs=DP('g',"level1"),
        theA=DP('a',"base")
  log:  **(EFS('allB.txt'))
  output:
    T("allB.txt")
  shell:
    "(time (echo 'the inputs: ' {input} > {output};        \
            echo 'allGs:' {params.allGs} >> {output}; \
            echo 'theA: ' {params.theA} >> {output} \
           )  2> {log.E}                    \
     ) 2> {log.T}"

rule level2_obj:
  input: 
    DT("obj.flag"),
    T('allB.txt')
  output:
    touch(T("obj.flag"))



