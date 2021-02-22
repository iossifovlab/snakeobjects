add_targets("allPairNumbers.txt","pairNumber.png")

rule gatherPairNumbers:
    input:  DT("pairNumber.txt")
    output: T("allPairNumbers.txt")
    shell:  "cat {input} | sort > {output}"

rule pairNumberFigure:
    input: T("allPairNumbers.txt")
    output: T("pairNumber.png")
    run:
        import matplotlib
        matplotlib.use("Agg")
        import pandas as pd 
        import matplotlib.pyplot as plt

        T = pd.read_table(input[0], sep='\t',header=None)

        fig = plt.figure(figsize=(3+len(T)/10, 3))
        plt.plot(T[1],'.')
        plt.xticks(range(len(T)),T[0],rotation=90.,fontsize=8)
        plt.xlim([-1,len(T)])
        plt.tight_layout()
        plt.savefig(output[0])
