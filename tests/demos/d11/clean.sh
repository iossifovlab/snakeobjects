#!/bin/bash
cd project_1
echo `pwd`
sobjects cleanProject -f
cd ../project_2
echo `pwd`
sobjects cleanProject -f
cd ../project_3
echo `pwd`
sobjects cleanProject -f
cd ../project_4
echo `pwd`
sobjects cleanProject -f

