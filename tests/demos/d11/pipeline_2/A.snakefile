add_targets("c.txt")
a = PP("P2_a")
b = PP("P2_b")
c = PP("P2_c")
d = PP("P2_d")

rule c:
  output: T("c.txt")
  shell: "echo {a} {b} {c} {d} >{output}"


rule e:
  output: T("e.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"

