add_targets("chrAll.bwaIndex.flag")

rule make_bwa_index:
    input: T("chrAll.fa")
    output:
         T("chrAll.bwaIndex.flag"),
    conda: "env-bwa.yaml"
    resources: mem_mb=10*1024
    log: **LFS("bwa_index")
    shell: "(time bwa index {input} -a bwtsw > {log.O} 2> {log.E} && touch {output[0]}) 2> {log.T}"
