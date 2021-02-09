rule fastqSummary:
  input:
    DT("obj.flag")
  output:
    touch(T("obj.flag"))

rule cnt:
  input:
    DT("sample_cnt.txt")
  output:
    T("sample_cnt.txt")
  log:
    **(LFS("sample_cnt.txt"))
  shell:
    "$SO_PIPELINE/summary_cnt.py {input} >{output} 2>{log.E}"

rule clean_fastqSummary:
  shell:
    " cd fastqSummary/o; rm -rf *"

