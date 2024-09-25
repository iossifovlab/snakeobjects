#!/bin/bash
set -e
cd $(dirname "$0")
sobjects prepare 5 
sobjects run -j 1
ls objects/P/*/obj.flag | wc

sobjects prepare 6 
sobjects run -j 1
ls objects/P/*/obj.flag | wc
