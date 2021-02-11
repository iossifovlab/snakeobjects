add_targets('allB.txt', 'merged.txt')

rule l2_allB:
  input:
    DT('B2.txt',"level1"),
    DT('a2.txt',"base")
  params: 
        allGs=DP('g',"level1"),
        theA=DP('a',"base")
  log:  **(LFS('allB.txt'))
  output:
    T("allB.txt")
  shell:
    "(time (echo 'the inputs: ' {input} > {output};        \
            echo 'allGs:' {params.allGs} >> {output}; \
            echo 'theA: ' {params.theA} >> {output} \
           )  2> {log.E}                    \
     ) 2> {log.T}"

chrs = ['1', '2', '3']

rule beg:
  output:
    T('beg-{c}.txt')
  shell:
    "echo 'this is beg' {wildcards.c} > {output}"

rule part:
  input:
    T('beg-{c}.txt')
  output:
    T('part-{c}.txt')
  shell:
    "echo 'this is part' {wildcards.c} > {output}; "
    "echo {input} >> {output}"

rule mergedI:
  input:
    expand(TE('part-{c}.txt'),c=chrs)
  output:
    T('merged.txt')
  shell:
    "cat {input} > {output}"


