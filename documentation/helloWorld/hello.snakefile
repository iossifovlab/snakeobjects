add_targets("result.txt")
rule createResult:
   output: T("result.txt")
   shell: "echo 'hello world' > {output}"
