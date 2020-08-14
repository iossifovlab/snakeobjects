rule P_big:
  input: 
    DT("a.txt")
  output:
    T('B.txt')
  log:  **(EFS('B.txt'))
  params: 
    thea = DP('a'),
    name = P('name'),
    dob  = P('dob'),
    ref  = config["parameters"]['ref']
  run:
    assert input[0] == "B/o/a.txt"
    correct = {
      "P/1/B.txt": ['Peter', "3/10/2000"],
      "P/2/B.txt": ['Paul' , "4/20/2001"],
      "P/3/B.txt": ['Mary' , "5/30/2002"],
      "P/4/B.txt": ['John' , "6/11/2003"]
    }
    assert output[0] in correct
    assert params.thea[0] == "alabala nica"
    assert params.name == correct[output[0]][0]
    assert params.dob  == correct[output[0]][1]
    shell("(time echo 'thea:' {params.thea}, 'name:' {params.name}, \
          'dob: {params.dob}'> {output}) 2> {log.T}")

rule P_obj:
  input: 
    DT("obj.flag"),
    T('B.txt')
  output:
    touch(T("obj.flag"))
  run:
    correct = {
    "P/1/obj.flag": ["B/o/obj.flag", "P/1/B.txt"],
    "P/2/obj.flag": ["B/o/obj.flag", "P/2/B.txt"],
    "P/3/obj.flag": ["B/o/obj.flag", "P/3/B.txt"],
    "P/4/obj.flag": ["B/o/obj.flag", "P/4/B.txt"]
    }
    assert output[0] in correct
    assert input == correct[output[0]]


