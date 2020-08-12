rule big_report:
  input:
    DT('realigned.bam','sample'),
    DT('ref.fa','index'),
  output:
    T('report.txt')
  params:
    allBams=DP('bamFile','sample')
  shell:
    "echo {input} > {output} ;"
    "echo {params.allBams} >> {output}"

rule report_obj:
  input:
    DT('obj.flag'),
    T('report.txt')
  output: 
    touch(T('obj.flag'))
