rule F_obj:
  input:
    DT('obj.flag')
  output: 
    touch(T('obj.flag'))
  params:
    state = P('state')
  run:
    correct = {
      "F/1/obj.flag": ["happy", "P/1/obj.flag", "P/2/obj.flag"],
      "F/2/obj.flag": ["  sad", "P/3/obj.flag", "P/4/obj.flag"]
    }
    assert output[0] in correct
    assert params.state == correct[output[0]][0]
    assert input[0] == correct[output[0]][1]
    assert input[1] == correct[output[0]][2]


