add_targets("sample_cnt.txt")

rule cnt:
  input:
    DT("sample_cnt.txt")
  output:
    T("sample_cnt.txt")
  log:
    **(LFS("sample_cnt.txt"))
  shell:
    "summary_cnt.py {input} >{output} 2>{log.E}"

rule clean_fastqSummary:
  shell:
    " cd fastqSummary/o; rm -rf *"

