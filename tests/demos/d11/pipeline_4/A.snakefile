add_targets("c.txt", "e.txt")
a = PP("P4_a")
b = PP("P4_b")
c = PP("P4_c")
d = PP("P4_d")

rule c:
  output: T("c.txt")
  shell: "echo {a} {b} {c} {d} >{output}"

rule e:
  output: T("e.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"

