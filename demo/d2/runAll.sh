for p in proj*; do
    (cd $p; 
        . ./setenv.sh
        pwd
        echo PROJECT_DIR $PROJECT_DIR
        echo PIEPELINE_DIR $PIPELINE_DIR
        build_object_graph.py createDirs
        run_snake.sh -j 
    )
done

