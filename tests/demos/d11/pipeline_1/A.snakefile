add_targets("a.txt","b.txt")
p = PP("project_1_A")
rule a:
  output: T("a.txt")
  shell: "a.sh {p}':' 'Once upon a time' >{output}"

rule b:
  output: T("b.txt")
  shell: "b.py 'there lived an old man and an old woman' >{output}"

