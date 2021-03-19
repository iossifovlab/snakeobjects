add_targets("sample.bam", "sample.bam.bai")
rule align:
  input:  ref=DT("ref.fa"), idx=DT("index.flag")
  output: temp(T("fastq.bam"))
  params: fD=PP('fastqDirectory'), fI=P('fqId')
  shell:  "bwa mem -R '@RG\\tID:{wildcards.oid}\\tSM:{wildcards.oid}' {input.ref}  \
                     {params.fD}/{params.fI}_1.fqz {params.fD}/{params.fI}_2.fqz | \
           samtools view -Sb - > {output}"
rule sort:
  input: T("fastq.bam")
  output: T("sample.bam")
  shell: "samtools sort {input} -O bam > {output}"
rule index:
  input: T("sample.bam")
  output: T("sample.bam.bai")
  shell: "samtools index -b {input} {output}"
