. setenv.sh
cp P1 P.snakefile
./build_object_graph.py createDirs
run_snake.sh -j -R 
ls -ltr objLinks/P/* >test_out.txt

echo "test dependency of obj.flag on new target" >>test_out.txt
cp P2 P.snakefile
run_snake.sh -j -R P_c 
ls -ltr objLinks/P/*  >>test_out.txt

./clean.sh
cp P1 P.snakefile
./build_object_graph.py createDirs
run_snake.sh -j -R 
ls -ltr objLinks/P/*  >>test_out.txt

echo "test dependency of old target on the new one"  >>test_out.txt
cp P3 P.snakefile
run_snake.sh -j -R P_c 
ls -ltr objLinks/P/*  >>test_out.txt


