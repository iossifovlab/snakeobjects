rule B_obj:
  input:
    DT('obj.flag')
  output: 
    touch(T('obj.flag'))
