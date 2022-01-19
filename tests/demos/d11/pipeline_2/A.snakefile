add_targets("a.txt", "b.txt", "c.txt")
a = PP("project_2_A")

rule a:
  output: T("a.txt")
  shell: "b.sh {a} >{output}"


rule b:
  output: T("b.txt")
  shell: " d.py `echo 'PYTHONPATH' $PYTHONPATH '\n\n'` >{output} && echo 'PATH' $PATH >>{output}"

rule c:
  output: T("c.txt")
  shell: "a.sh {a} >{output}"

