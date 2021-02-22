add_targets("sample.bam", \
            "sample_ns.bam", \
            "sample_cs.bam", "sample_cs.bam.bai", \
            "sample_cnt.txt")

rule bwa_map:
  input:
    ref=DT("chrAll.fa"),
    indexes=DT("obj.flag")
  output:
    T("sample.bam")
  conda:
    "envs/bwa.yaml"
  log:
    **(LFS('sample.bam'))
  params:
    r1 = P('R1'),
    r2 = P('R2'),
    rg = P('rg')
  shell:
    "bwa mem -R '{params.rg}' {input.ref} {params.r1} {params.r2} 2> {log.E} | samtools view -Sb - > {output}"

rule fastq_sort:
  input:
    T("sample.bam")
  output:
    T("sample_ns.bam")
  conda:
    "envs/bwa.yaml"
  shell:
        "samtools sort -n -T {input} -O bam {input} > {output}"

rule fastq_cs:
  input:
      T("sample_ns.bam")
  output:
    T("sample_cs.bam")
  log:
    **(LFS("sample_cs.bam"))
  conda:
    "envs/bwa.yaml"
  shell:
    "samtools sort -T {input} -O bam {input} > {output} 2>{log.E}"

rule fastq_idx:
  input:
    T("sample_cs.bam")
  output:
    T("sample_cs.bam.bai")
  conda:
    "envs/bwa.yaml"
  shell:
    "samtools index -b {input} {output}"

rule read_cnt:
  input:
    bam=T("sample_cs.bam"),
    idx=T("sample_cs.bam.bai")    
  output:
    T("sample_cnt.txt")
  params:
    regions=PP("target"),
    ref = PP("ref")
  conda:
    "envs/bwa.yaml"
  log:
    **(LFS("sample_cnt.txt"))
  shell:
    # -j option does not work for samtools 1.7
    "samtools bedcov --reference {params.ref} {params.regions} {input.bam} > {output} 2>{log.E}"

rule clean_fastq:
  shell:
    "cd fast; for n in *; do (cd $n; rm -rf *) done"
