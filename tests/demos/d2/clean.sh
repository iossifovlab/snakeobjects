#!/bin/bash
cd $(dirname "$0")
for p in proj*; do
    (cd $p; rm -rf *.snakefile OG.json OG.OG glbl.makefile log objects base)
done
