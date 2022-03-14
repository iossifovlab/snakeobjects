#!/bin/bash
cd $(dirname "$0")
rm -rf Snakefile && sobjects cleanProject -f
