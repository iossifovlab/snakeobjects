rule F_b:
  input:
    DT('b.txt', level=2, mode='lessOrEqual')
  output: 
    T('b.txt')
  params:
    DP('name', "P")
  run:
    correct = {
        "F/1/b.txt": ["P/1/b.txt", "P/2/b.txt", "B/o/b.txt", "Peter", "Paul"],
        "F/2/b.txt": ["P/3/b.txt", "P/4/b.txt", "B/o/b.txt", "Mary", "John"]
    }
    o = str(output)
    p1 = eval(str(params))[0]
    p2 = eval(str(params))[1]
    # print("OOOOOO",o,correct[o],input,p1, p2,file=sys.stderr)
    assert o in correct
    assert input == correct[o][:3]
    assert p1 == correct[o][3]
    assert p2 == correct[o][4]
    shell("echo {input} > {output}")


rule F_obj:
  input:
    DT('obj.flag'),
     T('b.txt')
  output: 
    touch(T('obj.flag'))
  params:
    state = P('state')
  run:
    correct = {
        "F/1/obj.flag": ["happy", "P/1/obj.flag", "P/2/obj.flag", "F/1/b.txt"],
        "F/2/obj.flag": ["  sad", "P/3/obj.flag", "P/4/obj.flag", "F/2/b.txt"]
    }
    assert output[0] in correct
    assert params.state == correct[output[0]][0]
    assert input[0] == correct[output[0]][1]
    assert input[1] == correct[output[0]][2]
    assert input[2] == correct[output[0]][3]


