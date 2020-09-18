localrules: split
chunkN = 7

rule P_b:
  input: 
    DT("a.txt")
  output:
    T('b.txt')
  log:  **(EFS('b.txt'))
  params: 
    thea = DP('a'),
    name = P('name'),
    dob  = P('dob'),
    ref  = GP('ref')
  run:
    assert input[0] == "B/o/a.txt"
    logs = ['log/b.txt-' + x for x in "out.txt err.txt time.txt".split(" ")]
    correct = {
        "P/1/b.txt": ['Peter', "3/10/2000"]+['P/1/' + x for x in logs],
        "P/2/b.txt": ['Paul' , "4/20/2001"]+['P/2/' + x for x in logs],
        "P/3/b.txt": ['Mary' , "5/30/2002"]+['P/3/' + x for x in logs],
        "P/4/b.txt": ['John' , "6/11/2003"]+['P/4/' + x for x in logs]
    }
    assert output[0] in correct
    assert params.thea[0] == "alabala nica"
    assert params.name == correct[output[0]][0]
    assert params.dob  == correct[output[0]][1]
    assert params.ref == "ref.fa"
    assert log.O == correct[output[0]][2]
    assert log.E == correct[output[0]][3]
    assert log.T == correct[output[0]][4]
    shell("(time echo 'thea:' {params.thea}, 'name:' {params.name}, \
          'dob: {params.dob}'> {output}) 2> {log.T}")


P_list = ["%05.i" % k for k in range(chunkN)]

rule split:
  input: 
    DT("a.txt", "B")
  output:
    touch(T("split.flag"))
  shell:
    " split -a 5 -d -n l/{chunkN} {input} `dirname {output}`/part- "

rule part:
  input:
    T("split.flag")
  output:
    T('part-{p}.txt')
  shell:
    "(echo 'this is part' {wildcards.oid} {wildcards.p}; cat `dirname {input}`/part-{wildcards.p}) > {output} && cp {output}  `dirname {input}`/c.txt "

rule mergedI:
  input:
    expand(TE('part-{p}.txt'),p=P_list)
  output:
    m=T('merged.txt'),
    c=touch(T("c.txt"))
  shell: 
    "cat {input} > {output.m} "
   
rule P_obj:
  input: 
    DT("obj.flag"),
    T('b.txt'),
    T('merged.txt')
  output:
    touch(T("obj.flag"))
  log:
    **(EFS("obj.flag"))
  run:
    correct = {
    "P/1/obj.flag": ["B/o/obj.flag", "P/1/b.txt", "P/1/merged.txt"],
    "P/2/obj.flag": ["B/o/obj.flag", "P/2/b.txt", "P/2/merged.txt"],
    "P/3/obj.flag": ["B/o/obj.flag", "P/3/b.txt", "P/3/merged.txt"],
    "P/4/obj.flag": ["B/o/obj.flag", "P/4/b.txt", "P/4/merged.txt"]
    }
    assert output[0] in correct
    assert input == correct[output[0]]
