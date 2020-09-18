rule AF_obj:
  input:
    DT('obj.flag')
  output: 
    touch(T('obj.flag'))
  resources: 
    mem_mb=500
  run:
    assert output[0] == "AF/o/obj.flag"
    for i in range(2):
      assert input[i] == "F/"+str(i+1)+"/obj.flag"
      assert input[2] == "B/o/obj.flag" 
