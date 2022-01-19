add_targets("a.txt","b.txt", "c.txt")
name = PP("project_1_A")

rule a:
  output: T("a.txt")
  shell: "a.sh $1 {name} >{output}"

rule b:
  output: T("b.txt")
  shell: "b.py `pwd`>{output}"

rule c:
  output: T("c.txt")
  shell: " echo 'PATH' $PATH '\n\n' >{output} && echo 'PYTHONPATH' $PYTHONPATH >>{output}"



