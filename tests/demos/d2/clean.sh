#!/bin/bash
cd $(dirname "$0")
rm pipeline/Snakefile
for p in proj*; do
    (cd $p; sobjects cleanProject -f)
done

