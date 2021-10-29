#!/bin/bash
set -e
cd $(dirname "$0")
cd project_1
sobjects prepare
sobjects run -j
cd ../project_2
sobjects prepare
sobjects run -j
