add_targets('assert.flag')

rule AP_assert:
  input:
    DT('obj.flag')
  output: 
    touch(T('assert.flag'))
  resources: 
    mem_mb=500
  run:
    assert output[0] == "AP/o/assert.flag"
    for i in range(4):
      assert input[i] == "P/"+str(i+1)+"/obj.flag"

