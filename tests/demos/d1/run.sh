#!/bin/bash
set -e
cd $(dirname "$0")
echo "HHHHHH" $0
echo "GGGGGG" `pwd`

echo 'RUNNING sobjects prepare', `type sobjects`
sobjects prepare

echo 'RUNNING sobjects run', `type sobjects`
sobjects run -j
