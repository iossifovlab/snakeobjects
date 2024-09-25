#!/bin/bash
set -e
cd $(dirname "$0")
for p in proj*; do
    (cd $p; 
        sobjects prepare
        sobjects run -j 1
    )
done

