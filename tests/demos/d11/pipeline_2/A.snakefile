add_targets("c.txt","d.txt")
p = PP("A")
a = PP("a")
b = PP("b")
rule c:
  output: T("c.txt")
  shell: "cat {a} >{output} && c.sh {p}':' 'I am here too' >>{output}"

rule d:
  output: T("d.txt")
  shell: "cat {b} > {output} && b.py 'coo coo' >>{output}"

rule e:
  output: T("e.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"

