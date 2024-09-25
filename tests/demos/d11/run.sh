#!/bin/bash
set -e
cd $(dirname "$0")
cd project_1
sobjects prepare
sobjects run -j 1
cd ../project_2
sobjects prepare
sobjects run -j 1
cd ../project_3
sobjects prepare
sobjects run -j 1
cd ../project_4
sobjects prepare
sobjects run -j 1
