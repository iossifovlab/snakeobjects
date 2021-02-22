add_targets('denovo_calls.vcf',"denovo_annotated.vcf")

rule summary_denovos:
  input:
    DT("denovo_calls.vcf")
  output:
    T("denovo_calls.vcf")
  shell:
    " head -5 {input[0]} >{output} && for n in {input}; do (tail -n +6 $n) done |sort -n -k3 >> {output} "

rule annotate_denovos:
  input:
    T("denovo_calls.vcf")
  output:
    T("denovo_annotated.vcf")
  conda:
    "envs/summary.yaml"
  log:
    **(LFS("denovo_annotated.vcf"))
  shell:
    " annotate.py {input} >{output}"
