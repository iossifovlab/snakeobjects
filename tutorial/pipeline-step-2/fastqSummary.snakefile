rule fastqSummary:
  input:
    T("allReadNumbers.txt"),
    T("readNumbmer.png"),
    DT("obj.flag")
  output:
    touch(T("obj.flag"))

rule gatherReadNumbers:
    input:  DT("readNumber.txt")
    output: T("allReadNumbers.txt")
    shell:  "cat {input} | sort > {output}"

rule readNumberFigure:
    input: T("allReadNumbers.txt")
    output: T("readNumbmer.png")
    run:
        import matplotlib
        matplotlib.use("Agg")
        import pandas as pd 
        import matplotlib.pyplot as plt

        T = pd.read_table(input[0], sep='\t',header=None)

        fig = plt.figure(figsize=(40, 3))
        plt.plot(T[1],'.')
        plt.xticks(range(len(T)),T[0],rotation=90.,fontsize=8)
        plt.xlim([-1,len(T)])
        plt.tight_layout()
        plt.savefig(output[0])
