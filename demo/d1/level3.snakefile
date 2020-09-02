rule level3:
  input:
    DT('obj.flag')
  output:
    touch(T('obj.flag'))

