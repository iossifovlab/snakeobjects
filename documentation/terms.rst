Glossary
---------

.. glossary::

    pipeline
        The pipeline usually contains a python script called ``build_object_graph.py`` that uses *meta data* associated with the projects that use the pipeline to create project's object graph and  a ``<object type>.snakefile`` for each of the object types created by the ``build_object_graph.py``. It also may contain executable files: python scripts, bash scripts, etc. as well as environment yaml files used in the workflow.

    project
        Directory containing ``so_project.yaml`` file with definition of parameters used in building object graph and in executing the pipeline.


    object graph
        `Object graph` is a structure representing a directed acyclic graph of *objects*  (the :py:class:`.ObjectGraph` is the ``snakeobjects`` implementation of the *object graph* and the objects in the object graph are implemented by the :py:class:`.OGO` class).

    objects
        A subdirectory of project directory. Its subdirectories are .snakeobjects, .snakemake, and a subdirectory for each object type in the object graph.
	
    target
        Target is a file associated with an object and resides in the object directory.

    object
        Object is an entity in the snakeobjects' object graph. It has associated directory, a list of objects it depends on, and a list of parameters.

    environment.yaml
        A file in the snakeobjects package specifying the minimal conda environment needed to successfully run snakeobjects projects.

    so_project.yaml
        A file stored in the project directory containing parameters that specify the pipeline operating on the project, pointers to the input and metadata associated with the project, and other parameters that control the processing.

    build_object_graph.py
        A python script that uses meta data associated with the projects to create project’s object graph and  `<object type>.snakefiles` templates for each of the object types that do not have already `<object type>.snakefile` in the pipeline directory. It also creates a directory `objects` in project directory and files `objects/.snakeobjects/main.snakefile`, `objects/.snakeobjects/OG.json` and directories for all objects in object graph.

    OG.json
        File containing json representation of `object graph`. It is created by command :option:`sobjects prepare`.

    main.snakefile
        File passed to snakemake by command :option:`sobjects run`. It is created by command :option:`sobjects prepare`.

    jobscript.sh
        Bash script created by command :option:`sobjects submit` and passed to corresponding cluster engine command (i.e., 'sbash' for slurm or 'qsub' for sge). Its location is in `objects/.snakeobjects`.

    .snakeobjects
        Subdirectory of `objects` directory containing files `jobscript.sh`, `OG.json` and `main.snakefile`.

    .snakemake
        Subdirectory of `objects` directory that is created by running snakemake. 

    sobjects commands
        See full list of command in :ref:`sobjects-commands`.

    object type
        An `object type` is characterized by a set of `targets` that need to be created for each `object` of the given object type together with the rules for creating the targets.

    object parameters
        Parameters derived from project `metadata`, specific to concrete object and accessible in snakemake rules.

    project parameters
        Parameters derived from project `metadata` common to all objects and accessible in snakemake rules.

    object name
        Each object has a proper name. A name is a character string satisfying conventions of file names and json entities names. 

    meta data, or metadata
        The project `meta data` can have an arbitrary form (a list input files; a csv file; relation database, etc.) and is usually used to generate the project-specific object graph.

