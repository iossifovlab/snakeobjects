add_targets('assert.flag')

rule AF_assert:
  input:
    DT('obj.flag')
  output: 
    touch(T('assert.flag'))
  resources: 
    mem_mb=500
  run:
    assert output[0] == "AF/o/assert.flag"
    for i in range(2):
      assert input[i] == "F/"+str(i+1)+"/obj.flag"
      assert input[2] == "B/o/obj.flag" 
