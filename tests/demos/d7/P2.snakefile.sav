add_targets("a.txt","b.txt","c.txt")

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
  output: 
    touch(T('b.txt'))

