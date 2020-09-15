. ./setenv.sh
./build_object_graph.py 5 createDirs
run_snake.sh -j
ls objLinks/P/*/obj.flag | wc
./build_object_graph.py 6 createDirs
run_snake.sh -j
ls objLinks/P/*/obj.flag | wc
