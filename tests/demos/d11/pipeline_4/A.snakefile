add_targets("c.txt", "e.txt")
a = PP("P4_a")
b = PP("P4_b")


rule c:
  output: T("c.txt")
  shell: "c.sh {a} {b} >{output}"

rule e:
  output: T("e.txt")
  shell: " d.py `echo 'PATH' $PATH '\n\n'` >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"

