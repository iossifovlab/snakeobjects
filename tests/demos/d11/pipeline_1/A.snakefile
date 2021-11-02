add_targets("a.txt","b.txt", "e.txt")
a = PP("P1_a")
b = PP("P1_b")

rule a_:
  output: T("a.txt")
  shell: "echo {a} >{output}"

rule b_:
  output: T("b.txt")
  shell: "echo {b} >{output}"

rule e:
  output: T("e.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"


