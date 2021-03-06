add_targets('assert.flag')

rule P_b:
  input: 
    DT("a.txt")
  output:
    T('b.txt')
  log:  **(LFS('b.txt'))
  params: 
    thea = DP('a'),
    name = P('name'),
    dob  = P('dob'),
    ref  = PP('ref')
  resources: 
    mem_mb=500
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


chrs = ['chr1', 'chr2', 'chr3', 'chr4']

rule beg:
  output:
    T('beg-{c}.txt')
  resources: 
    mem_mb=500
  shell:
    "echo 'this is beg' {wildcards.oid} {wildcards.c} > {output}"

rule part:
  input:
    T('beg-{c}.txt')
  output:
    T('part-{c}.txt')
  resources: 
    mem_mb=500
  shell:
    "echo 'this is part' {wildcards.oid} {wildcards.c} > {output}; "
    "echo {input} >> {output}"

rule mergedI:
  input:
    expand(TE('part-{c}.txt'),c=chrs)
  output:
    T('merged.txt')
  resources: 
    mem_mb=500
  shell: 
    "cat {input} > {output}"


rule P_assert:
  input: 
    DT("obj.flag"),
    T('b.txt'),
    T('merged.txt')
  output:
    touch(T("assert.flag"))
  resources: 
    mem_mb=500
  run:
    correct = {
    "P/1/assert.flag": ["B/o/obj.flag", "P/1/b.txt", "P/1/merged.txt"],
    "P/2/assert.flag": ["B/o/obj.flag", "P/2/b.txt", "P/2/merged.txt"],
    "P/3/assert.flag": ["B/o/obj.flag", "P/3/b.txt", "P/3/merged.txt"],
    "P/4/assert.flag": ["B/o/obj.flag", "P/4/b.txt", "P/4/merged.txt"]
    }
    assert output[0] in correct
    assert input == correct[output[0]]


