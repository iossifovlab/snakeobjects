rule trioSummary:
  input:
    DT("obj.flag"),
    T("denovo_calls.txt")
  output:
    touch(T("obj.flag"))

rule summary_denovos:
  input:
    DT("denovo_calls.txt")
  output:
    T("denovo_calls.txt")
  shell:
    " head -1 {input[0]} >{output} && for n in {input}; do (tail -n +2 $n) done |sort -n -k3 >> {output} "

