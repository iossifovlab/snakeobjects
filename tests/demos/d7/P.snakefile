rule P_obj:
  input:
    DT('obj.flag'),
     T('a.txt'),
     T('b.txt')
  output: 
    touch(T('obj.flag'))

rule P_a:
  input:
    DT('obj.flag')
  output: 
    touch(T('a.txt'))

rule P_c:
  input:
    T('a.txt')
  output: 
    touch(T('c.txt'))
	 
rule P_b:
  input:
    T('a.txt'),
    T('c.txt')
  output: 
    touch(T('b.txt'))

