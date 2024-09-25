#!/bin/bash
set -e
cd $(dirname "$0")
sobjects prepare
sobjects run -j 1 --use-conda --conda-frontend conda 
