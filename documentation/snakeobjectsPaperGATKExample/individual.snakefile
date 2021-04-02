add_targets('sample.bam', 'sample.bam.bai', 'output.g.vcf')
rule align:
  input:  ref=DT('ref.fa'), idx=DT('index.flag')
  output: temp(T('fastq.bam'))
  params: fD=PP('fastqDirectory'), fI=P('fqId')
  shell:  "bwa mem -R '@RG\\tID:RG1\\tSM:{wildcards.oid}' \
              {input.ref} {params.fD}/{params.fI}_1.fqz   \
                          {params.fD}/{params.fI}_2.fqz | \
           samtools view -Sb - > {output}"
rule sort:
  input:  T('fastq.bam')
  output: T('sample.bam')
  shell:  "samtools sort {input} -O bam > {output}"

rule index:
  input:  T('sample.bam')
  output: T('sample.bam.bai')
  shell:  "samtools index -b {input} {output}"

rule haplotypeCaller:
  input:  Bs=T('sample.bam'), Is=T('sample.bam.bai'), ref=DT('ref.fa') 
  output: T('output.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' HaplotypeCaller -ERC GVCF -R {input.ref} -I {input.Bs} -O  {output}"

rule genotypeCaller:
  input:  vcf=T('output.g.vcf'), ref=DT('ref.fa') 
  output: T('genotype.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' GenotypeGVCFs -R {input.ref} -V {input.vcf} -O  {output}"