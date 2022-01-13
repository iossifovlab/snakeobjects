add_targets("c.txt","e.txt")
a = PP("P3_a")
b = PP("P3_b")

rule c:
  output: T("c.txt")
  shell: "echo {a} {b} >{output}"

rule e:
  output: T("e.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"

