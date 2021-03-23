.. _getting-started-label:

***************
Getting Started
***************

Installation
============

conda install -c iossifovlab -c bioconda -c conda-forge snakeobjects

Minimal (empty) project and pipeline
====================================

Step 1
------

An empty directory is both a valid ``snakeobjects`` project and a valid
``snakeobjects`` pipeline directory.  We will use here such a minimal and
useless project to introduce some of the basic steps for working with
``snakeobjects`` projects and pipelines. Let's create an empty directory and go
in it:

.. code-block:: bash

    $ mkdir /tmp/minimalSO
    $ cd /tmp/minimalSO

Step 2
------

We can then use the :option:`sobjects describe` tool to get basic information about
the project and the pipeline:

.. code-block:: bash

    $ sobjects describe
    WORKING ON PROJECT /tmp/minimalSO
    WITH PIPELINE /tmp/minimalSO
    Project parameters:
    Types:

The result shows that we are using a project in the directory
``/tmp/minimalSO`` and that this project uses a pipeline in the same
``/tmp/minimalSO`` directory. Following, is an empty list of the *Project
parameters* (no ``so_project.yaml`` file is provided for the project and thus
there are no project parameters) and an empty list of *Object types* (no
objects have been added to the object graph so there are no object types used).

Step 3
^^^^^^

We then use the :option:`sobjects prepareObjects` command to prepare the projects for execution: 

.. code-block:: bash

    $ sobjects prepareObjects
    WORKING ON PROJECT /tmp/minimalSO
    WITH PIPELINE /tmp/minimalSO
    $ find .
    .
    ./objects
    ./objects/.snakeobjects
    ./objects/.snakeobjects/main.snakefile

This command assembles the snakefile to be used by ``snakemake`` and stores it
in the ``./objects/.snakeobjects/main.snakefile``.  The command would also
create directories for all objects in the object graph, but, since no objects
are added to the object graph, no object directories are created.

Step 4
^^^^^^

Finally, we can run :option:`sobjects run` command to create all targets for all the
object in the project:

.. code-block:: bash

    $ sobjects run -j -q
    WORKING ON PROJECT /tmp/minimalSO
    WITH PIPELINE /tmp/minimalSO
    RUNNING: snakemake -s /tmp/minimalSO/objects/.snakeobjects/main.snakefile -d /tmp/minimalSO/objects -j -q
    Job counts:
        count	jobs
        1	all_main
        1

The command shows the project and pipeline directories it will use and the
command line used to run ``snakemake``.  All parameters given to ``sobjects
run`` (i.e. ``-j -q``) are passed directly to ``snakemake``). ``-j`` instructs
``snakemake`` to run the pipeline on the local host and to use all the
available processors and  ``-q`` instruct ``snakemake`` to be *quiet*
(``snakemake`` is rather verbose by default). Since no objects are present in the object
graph, no targets are created. The only change is that ``snakemake`` creates
its own private directory in ``./objects/.snakemake``.

``hello world`` project
=======================

Here will show how to create a project with only one object with object type ``hello`` and object id ``world``. 
As above we will use the same directory for the project and for the pipeline:  

.. code-block:: bash

    $ mkdir /tmp/minimalHW
    $ cd /tmp/minimalHW

But now will add one object to the project's object graph. To do that we have to create a python 
file called ``build_object_graph.py`` 
in the pipeline directory with the following content:

.. code-block::

    def run(proj,OG):
        OG.add("hello","world")


