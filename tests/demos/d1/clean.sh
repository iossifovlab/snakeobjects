#!/bin/bash
cd $(dirname "$0")
rm -r Snakefile && sobjects cleanProject -f
