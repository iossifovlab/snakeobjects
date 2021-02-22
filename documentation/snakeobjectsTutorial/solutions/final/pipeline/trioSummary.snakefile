add_targets("allDenovoCalls.txt")

rule gatherDenovos:
    input: DT("denovo_calls.txt")
    output: T("allDenovoCalls.txt")
    shell: '''
        head -1 {input[0]} > {output}
        for t in {input}; do
            tail -n +2 $t >> {output}
        done
        '''
