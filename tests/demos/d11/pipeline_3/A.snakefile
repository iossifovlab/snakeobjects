add_targets("c.txt","e.txt")
a = PP("P3_a")

rule c:
  output: T("c.txt")
  shell: "echo {a} >{output}"

rule e:
  output: T("e.txt")
  shell: " c.sh 'PATH' $PATH '\n\n' >{output}"

