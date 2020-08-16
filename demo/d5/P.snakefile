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
    ref  = PP('ref')
  run:
    assert input[0] == "B/o/a.txt"
    correct = {
      "P/1/b.txt": ['Peter', "3/10/2000"],
      "P/2/b.txt": ['Paul' , "4/20/2001"],
      "P/3/b.txt": ['Mary' , "5/30/2002"],
      "P/4/b.txt": ['John' , "6/11/2003"]
    }
    assert output[0] in correct
    assert params.thea[0] == "alabala nica"
    assert params.name == correct[output[0]][0]
    assert params.dob  == correct[output[0]][1]
    assert params.ref == "ref.fa"
    shell("(time echo 'thea:' {params.thea}, 'name:' {params.name}, \
          'dob: {params.dob}'> {output}) 2> {log.T}")

rule P_obj:
  input: 
    DT("obj.flag"),
    T('b.txt')
  output:
    touch(T("obj.flag"))
  run:
    correct = {
    "P/1/obj.flag": ["B/o/obj.flag", "P/1/b.txt"],
    "P/2/obj.flag": ["B/o/obj.flag", "P/2/b.txt"],
    "P/3/obj.flag": ["B/o/obj.flag", "P/3/b.txt"],
    "P/4/obj.flag": ["B/o/obj.flag", "P/4/b.txt"]
    }
    assert output[0] in correct
    assert input == correct[output[0]]


