
rule level1_obj:
  output:
    T("obj.flag")
  shell:
    "touch {output}"



