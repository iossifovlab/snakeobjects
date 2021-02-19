add_targets("chrAll.bwaIndex.flag")

rule make_bwa_index:
    input: T("chrAll.fa")
    output:
         T("chrAll.bwaIndex.flag"),
         T("chrAll.fa.amb"),
         T("chrAll.fa.ann"),
         T("chrAll.fa.bwt"),
         T("chrAll.fa.pac"),
         T("chrAll.fa.sa")
    conda: "env-bwa.yaml"
    shell: "bwa index {input} -a bwtsw && touch {output[0]}"
