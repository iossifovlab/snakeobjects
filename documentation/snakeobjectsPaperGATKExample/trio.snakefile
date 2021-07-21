add_targets('output.g.vcf', 'genotype.g.vcf')

rule combineGVCFs:
  input:  vcf=DT('output.g.vcf'), ref=DT('ref.fa', level=2) 
  output: T('output.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' CombineGVCFs -R {input.ref} `for x in {input.vcf}; do echo -V $x; done` -O  {output}"


rule genotypeGVCFs:
  input:  vcf=T('output.g.vcf'), ref=DT('ref.fa', level=2) 
  output: T('genotype.g.vcf')
  resources: mem_mb=6000
  shell:  "gatk --java-options '-Xmx5G' GenotypeGVCFs -R {input.ref} -V {input.vcf} -O  {output}"
  