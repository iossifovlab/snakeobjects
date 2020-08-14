
rule build_a:
  output:
     t = T("t.txt")
  log:  **(EFS('t.txt'))
  params: a=P('a')
  run:
    assert output[0] == "B/o/t.txt"
    assert params.a == "alabala nica"
    shell("echo {params.a} > {output.t}")

rule build_ab:
  output:
    T("a.txt"), T("b.txt")
  log:  **(EFS('b.txt'))
  run:
    assert output[0] == "B/o/a.txt"
    assert output[1] == "B/o/b.txt"
    shell("touch {output[0]}")
    shell("touch {output[1]}")

rule base_obj:
  input:
    T("a.txt"), T("b.txt"), T("t.txt")
  output:
    touch(T("obj.flag"))
  run:
    assert input[0]  == "B/o/a.txt"
    assert input[1]  == "B/o/b.txt"
    assert input[2]  == "B/o/t.txt"
    assert output[0] == "B/o/obj.flag"
