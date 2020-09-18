rule AP_obj:
  input:
    DT('obj.flag')
  output: 
    touch(T('obj.flag'))
  resources: 
    mem_mb=500
  run:
    assert output[0] == "AP/o/obj.flag"
    for i in range(4):
      assert input[i] == "P/"+str(i+1)+"/obj.flag"

