add_targets("ref.fa","index.flag")
rule indexRef:
  input:  PP("reference")
  output: T("ref.fa"),touch(T("index.flag"))
  shell: "ln -s {input} {output[0]} && \
          bwa index {output[0]} -a bwtsw"
