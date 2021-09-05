#!/bin/bash
set -e
cd $(dirname "$0")

rm -f test_out.txt

echo "BASELINE" >>test_out.txt
./clean.sh
cp P1.snakefile.sav P.snakefile
sobjects prepare
sobjects run -j -R 
ls P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt


echo "test dependency of obj.flag on new target" >>test_out.txt
cp P2.snakefile.sav P.snakefile
sobjects run -j -R P_c 
ls P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt


echo "BASELINE" >>test_out.txt
./clean.sh
cp P1.snakefile.sav P.snakefile
sobjects prepare
sobjects run -j -R 
ls P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt

echo "test dependency of old target on the new one"  >>test_out.txt
cp P3.snakefile.sav P.snakefile
sobjects run -j -R P_c 
ls P/*  | grep -v log >>test_out.txt

diff test_out.txt correct_test_out.txt
