rule AF_obj:
  input:
    DT('obj.flag'),
    T("t.txt")
  output: 
    touch(T('obj.flag'))
  run:
    assert output[0] == "AF/o/obj.flag"
    for i in range(2):
      assert input[i] == "F/"+str(i+1)+"/obj.flag"
      assert input[2] == "B/o/obj.flag" 


rule create_t:
    output: T("t.txt")
    shell: """
           echo PATH is $PATH > {output}
           echo PYTHONPATH is $PYTHONPATH >> {output}
           echo SO_PROEJCT is $SO_PROJECT >> {output}
           echo SO_PIPELINE is $SO_PIPELINE >> {output}
           s1.py >> {output}
           """

