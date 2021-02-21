add_targets("sample.bam", "sample.bam.bai","markDupStats.txt")

rule merge:
    input: DT("fastq.bam")
    output: temp(T("raw.bam"))
    conda: "env-bwa.yaml"
    log: **LFS('merge')
    shell: "(time samtools merge -n {output} {input} > {log.O} 2> {log.E} ) 2> {log.T}"

rule reorganized_bam:
    input: T("raw.bam"), DT("chrAll.fa",level=2)
    output: T("sample.bam"), T("markDupStats.txt")
    conda: "env-bwa.yaml"
    log: **LFS('reorganize')
    shell: '''
        (time 
            samtools fixmate -m --reference {input[1]} -O bam {input[0]} - |
            samtools sort -T {input[0]} -O bam | 
            samtools markdup -T {input[0]} -O bam -s --reference {input[1]} - {output[0]} 2> {output[1]}
        ) 2> {log.T}
    '''

rule sample_idx:
    input: T("sample.bam")
    output: T("sample.bam.bai")
    conda: "env-bwa.yaml"
    log: **LFS('merge')
    shell: "(time samtools index -b {input} {output} > {log.O} 2> {log.E} ) 2> {log.T}" 

add_targets("depth.txt","coverage-stats.txt","coverage.png")

rule target_depth:
    input:
        bam=T("sample.bam"),
        idx=T("sample.bam.bai")
    output: T("depth.txt")
    params: target=PP("target")
    conda: "env-bwa.yaml"
    log: **LFS('target_depth')
    shell: "(time samtools depth -b {params.target} -q 30 -Q 30 {input.bam} > {output} 2> {log.E}) 2> {log.T}"

rule coverage_stats:
    input: T("depth.txt")
    output: T("coverage-stats.txt")
    run:
        import pandas as pd 
        D = pd.read_table(input[0], sep='\t',header=None)[2]
        coverageStr = ['%.1f' % (100*sum(D[D >= k])/sum(D)) for k in [1,10,20,40]]
        with open(output[0],"w") as OF:
            OF.write(f'{wildcards.oid}\t' + "\t".join(coverageStr) + "\n")

rule coverage:
    input: T("depth.txt")
    output: T("coverage.png")
    run:
        import matplotlib
        matplotlib.use("Agg")
        import pandas as pd 
        import matplotlib.pyplot as plt

        D = pd.read_table(input[0], sep='\t',header=None)[2]
        coverage = [100*sum(D[D >= k])/sum(D) for k in range(41)]
        fig = plt.figure(figsize=(5, 3))
        plt.plot(coverage,'b.-')
        plt.grid(1)
        xtcks = range(0,41,5)
        plt.xticks(xtcks, [f'{d}x' for d in xtcks])
        plt.xlim([0,40])
        ytcks = range(0,101,20)
        plt.yticks(ytcks, [f'{p}%' for p in ytcks])
        plt.ylim([0,100])
        plt.ylabel('percent of the target covered at >= X')
        plt.title(f"Coverage for sample {wildcards.oid}")
        plt.tight_layout()
        plt.savefig(output[0])
