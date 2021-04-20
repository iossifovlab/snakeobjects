add_targets('sample.bam', 'sample.bam.bai', 'output.g.vcf', 'sample_mdup.bam')
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

rule markDuplicates:
  input:  Bs=T('sample.bam'), Bi=T('sample.bam.bai')
  output: b=T('sample_mdup.bam'), m=T('marked_dup_metrics.txt')
  resources: mem_mb=6000
  shell:  "gatk MarkDuplicates \
  	  -I {input.Bs} \
	  -O {output.b} \
	  -M {output.m}"

rule index_mdub:
  input:  T('sample_mdup.bam')
  output: T('sample_mdup.bam.bai')
  shell:  "samtools index -b {input} {output}"

rule haplotypeCaller:
  input:  Bs=T('sample_mdup.bam'), i=T('sample_mdup.bam.bai'), ref=DT('ref.fa') 
  output: T('output.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' HaplotypeCaller \
  	  -ERC GVCF      \
	  -R {input.ref} \
	  -I {input.Bs}  \
	  -O  {output}"

rule genotypeCaller:
  input:  vcf=T('output.g.vcf'), ref=DT('ref.fa') 
  output: T('genotype.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' GenotypeGVCFs -R {input.ref} -V {input.vcf} -O  {output}"