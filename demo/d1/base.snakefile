
rule build_a:
  output:
    t = T("t.txt")
  shell:
    "touch {output.t}; "

rule build_ab:
  output:
    T("a.txt"), T("b.txt")
  shell:
    "touch {output[0]}; "
    "touch {output[1]}; "

rule base_obj:
  input:
    T("a.txt"), T("b.txt"), T("t.txt")
  output:
    T("obj.flag")
  shell:
    "touch {output}"



