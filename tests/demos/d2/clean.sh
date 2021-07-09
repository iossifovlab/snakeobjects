#!/bin/bash
cd $(dirname "$0")
rm pipeline/Snakefile
for p in proj*; do
    (cd $p; rm -rf .snakemake OG.json Snakefile log objects base level1 level2)
done

