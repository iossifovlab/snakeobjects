add_targets("allCoverageStats.txt","allCoverages.png")

rule gatherCoverageStats:
    input:  DT("coverage-stats.txt")
    output: T("allCoverageStats.txt")
    shell:  "cat {input} | sort -k4n > {output}"

rule allCoveragesFigure:
    input: T("allCoverageStats.txt")
    output: T("allCoverages.png")
    run:
        import matplotlib
        matplotlib.use("Agg")
        import pandas as pd 
        import matplotlib.pyplot as plt

        T = pd.read_table(input[0], sep='\t',header=None)

        fig = plt.figure(figsize=(3+len(T)/10, 3))
        plt.plot(T[1],'.',label='X >= 1')
        plt.plot(T[2],'.',label='X >= 10')
        plt.plot(T[3],'.',label='X >= 20')
        plt.plot(T[4],'.',label='X >= 40')

        plt.xticks(range(len(T)),T[0],rotation=90.,fontsize=8)
        plt.xlim([-1,len(T)])
        plt.ylim([0,100])
        plt.ylabel('percent of target covered at')
        plt.legend(loc='lower left')
        plt.tight_layout()
        plt.savefig(output[0])
