#!/bin/bash
cd $(dirname "$0")
. setenv.sh

rm test_out.txt

echo "BASELINE" >>test_out.txt
./clean.sh
cp P1.snakefile P.snakefile
./build_object_graph.py createDirs
run_snake.sh -j -R 
ls objLinks/P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt


echo "test dependency of obj.flag on new target" >>test_out.txt
cp P2.snakefile P.snakefile
run_snake.sh -j -R P_c 
ls objLinks/P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt


echo "BASELINE" >>test_out.txt
./clean.sh
cp P1.snakefile P.snakefile
./build_object_graph.py createDirs
run_snake.sh -j -R 
ls objLinks/P/*  | grep -v log >>test_out.txt
echo  >> test_out.txt
echo  >> test_out.txt

echo "test dependency of old target on the new one"  >>test_out.txt
cp P3.snakefile P.snakefile
run_snake.sh -j -R P_c 
ls objLinks/P/*  | grep -v log >>test_out.txt

diff test_out.txt correct_test_out.txt
