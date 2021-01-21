Introduction
++++++++++++


Overview
--------

**snakeobjects** is a workflow management framework based on ``snakemake`` that
uses an object-oriented abstraction of workflows. ``snakeobjects`` workflows
are easier to develop, to maintain and to adopt compared to the equivalent
workflows written in ``snakemake``, but inherit all the powerful features of
``snakemake``. These include the portability, efficient resource usage, the
large expressive power due to the tight python integration, and the large
community of the ``snakemake`` users. 

Workflow object-oriented abstraction
------------------------------------

The ``snakeobjects`` introduces an abstraction of workflows inspired by
object-oriented design that replaces the low-level input-output relationships
between files that are at the core of snakemakes rules. A *pipeline* (workflow)
in snakeobjects operates on *projects* and projects are composed of *objects*. The
objects within a project are connected with *dependency relationships* organized
in a directed acyclic graph called *object graph* in which each object has a list
of the objects it depends on (*dependency objects*). Each object has also a
specified *object type* and an object type is characterized by a set of *targets*
that need to be created for each object of the given object type together with
the rules for creating the targets. The rules for building targets for an
object type are included in a snakefile named after the object type and are written
using snakemake's syntax where the inputs and outputs specify targets instead of
files. Crucially, inputs can refer to targets in the current object and to
targets in objects the current object depends on as specified in the object
graph. Finally, projects and objects can be associated with a set of key-value
parameters.

Pipelines 
---------

In ``snakeobjects``, pipelines reside in a *pipeline directory*. The pipeline
directory and its content are created by the *workflow designer* and define the
workflow. The pipeline usually contains a python script called
``build_object_graph.py`` that uses *meta data* associated with the projects
that use the pipeline to create project's object graph and  a ``<object type>.snakefile`` for each of the
object types created by the ``build_object_graph.py``. 

``build_object_graph.py`` script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``build_object_graph.py`` script that should contain a function with the following 
interface:

.. py:function:: run(project, OG [,*args])

        Creates an object graph for the **project**.
    
        :param project: the snakeobjects project
        :type project: snakeobjects.Project 

        :param OG: the newly created object graph
        :type OG: snakeobjects.ObjectGraph
 
        :param \*args: command line arguments passed to :option:`sobjects prepare` or 
                       :option:`sobjects prepareTest`.
        :type \*args: list[str]


The :py:func:`run` function usually obtains the location of the project meta data through
``project.parameters``, loads the meta data, and uses it to create the
corresponding object graph using the :py:meth:`~snakeobjects.ObjectGraph.add`
method to add object to the ``OG``.

Object-type snakefiles
^^^^^^^^^^^^^^^^^^^^^^

Each object-type snakefiles contains
the list the targets for the object type and the rules for creating the targets.
The list of targets is defined with a rule for the creating *special* target ``T("obj.flag")``.
The ``T("obj.flag")`` target for an object will be created only after all the object type targets have been created.

For example, the following rule in the object type snakefile ``sample.snakefile``

.. code-block::

    rule sample_obj:
       input:
            T("A.txt"),
            T("B.txt"),
            DT("obj.flag")
       output: 
            touch(T("obj.flag"))

, decleares that each of the objects of type ``sample`` need two targets
created (``T("A.txt")`` and ``T("B.txt")``).  Also, the speclial target
``T("obj.flag")`` will be created only after, the ``T("obj.flag")`` targets are
created for all objects the objects the current object depends on. This is
specified by the inclusion of the ``DT("obj.flag")`` list of input targets.

All the rules, including the one for the ``T("obj.flag")`` special target,
are written using the ``snakemake``'s syntax (https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) 
and use a set of ``snakeobjects`` extension functions (see
:ref:`snake-extenssions`)  for referring to targets in the current object
(:py:func:`.T`); to targets in dependency objects (:py:func:`.DT`); to
parameters of the current object (:py:func:`.P`); to parameters of the
dependency objects (:py:func:`.DP`); and to global object parameters
(:py:func:`.PP`).  In addition, the pipeline directory can contain scripts,
conda environment definitions or other artefacts used by the pipeline. 

The example below demonstrates the main features of the ``snakeobjects`` rules:

.. code-block:: python

    rule create_B:
        input: a=T("A.txt"), r=DT("chrAll.fa",dot="reference")
        output: T("B.txt")
        parameters: g=P("gender")
        log: **LFS("B")
        shell: "some_command.py {input.a} {param.g} {input.ref} > {output} 2> {log.E}"

.. TODO: Add description of the example above.

The first two lines in this rule use functions :py:func:`.T`  and :py:func:`.DT`.

Parameters and log are definen by functions :py:func:`.P` and :py:func:`.LFS`.

Projects
--------

A project in ``snakeobjects`` is created by a *workflow user* to apply one
``snakeobjects`` pipeline.  A project is associated with a *project directory*
that usually contains a ``so_project.yaml`` file written by the *workflow user*
to configure the project.  The *workflow user* uses the ``sobjects``
command line tool to initialize (usually using the :option:`sobjects prepare`
command) and to execute (:option:`sobjects run`) the associated
*pipeline*.  The :option:`sobjects prepare`:
 
1. creates and object graph (using 
   the ``build_object_graph.py`` script from the *pipeline*)
   and stores it in the ``snakeobjects``'s private subdirectory ``.snakeobjects`` of the *project directory* 
   (``<projet directory>/.snakeobjects/OG.json``); 
2. creates an *object diretory* directory for
   each of the objects in the *object graph* in the ``objects`` subdirectory in the *project directory*; 
3. creates the ``<project directory>/.snakeobjects/main.snakefile`` that
   is subsequently used by ``snakemake``; and 
4. creates the symbolic links based for all object that have ``symlink.<name>`` parameters. 

The targets and the log files created during the execution of the pipeline (:option:`sobjects run`) are 
stored in the *object directories* in the ``objects`` subdirectory. 
In addition, ``snakemake`` creates it's own standard internal 
subdirectory ``.snakemake`` as a subdirectory the *project directory*.

``so_project.yaml`` file
^^^^^^^^^^^^^^^^^^^^^^^^

The ``so_project.yaml`` contains the *project parameters* that configure the
project and may include:

* a ``so_pipeline`` parameter that points to the *pipeline directory* for the
  pipeline that will operate on the project (a relative paths are relative 
  based on the project directory);
* parameters pointing to the input that will be used by the project; 
* parameters pointing to the meta-data describing the projects input; 
* a ``default_snakemake_params`` parameter that specifies the command line 
  arguments that are passed to ``snakemake`` at every invocation of 
  :option:`sobjects run`; 


Parameter values may contain expressions ``[E:<env_variable_name>]``, ``[C:<parameter>]``,
or ``[P:<project property>]``.
In the first case the expression is replaced by the value of environment variable called 
``env_variable_name``; in the second case the expression is replaced with the
value of parameter called ``parameter`` in the ``so_project.yaml`` file; in the
third case the expression is replaced with the project directory if ``project
property`` is ``projectDir.
  

``objects`` subdirectory
^^^^^^^^^^^^^^^^^^^^^^^^

The files related to ``snakeobjects`` targets have the following general name::

    <project directory>/objects/<object type>/<object id>/<target name>

For example, the target ``T("A.txt")`` of object of object type ``sample`` and with
id ``i1232`` will be stored in the files ``<project
directory>/objects/sample/i1232/A.txt``; 


The general form for the ``log.O``, ``log.E``, and ``log.T`` log files referenced 
using the ``logEFS(<name>)`` function are::

    <project directory>/objects/<object type>/<object id>/log/<name>-out.txt
    <project directory>/objects/<object type>/<object id>/log/<name>-err.txt
    <project directory>/objects/<object type>/<object id>/log/<name>-time.txt

respectively. For example, log file (``log.E``)
named ``A`` for the sample i1232 object is ``<project
directory>/objects/sample/i1232/log/A-err.txt``. 

``.snakeobjects`` subdirectory
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is a private directory reserved for ``snakeobjects`` internal files. Currently, 
the directory contain two files that may be of interest to the *workflow user*:

* ``.snakeobjets/OG.json`` contains the object graph associated with the project;
* ``.snakeobjets/main.snakefile`` contains the snakefile that is passed to ``snakemake`` at the 
  :option:`sobjects run`. 

Objects types, objects, and object graph
----------------------------------------

Object types in ``snakeobjects`` typically correspond to object types in the domain of the 
*workflow*. For example, in sequence analysis *workflows*, we can have *reference genome*, 
*library*, *sample*, *individual*, *family*, or *population* object types. ``snakeobjects`` object
types are characterized by the set of *targets* that will be created for each object of the 
object type. For example, *sample* may have a ``T("sample.bam")``, 
a ``T("sample.bai")``, a ``T("sample.vcf")``, and a ``T("depth-histogram.png")``; 
*reference genome* object type may have a ``T("chr.fa")``, ``T("chr.fa.fai")``, 
``T("bwa.index")`` targets. 

Each ``snakeobjects`` project is associated with one :term:`object graph`
a structure representing 
a directed acyclic graph of *objects*  (the :py:class:`.ObjectGraph` is the ``snakeobjects`` implementation of the *object graph* and the objects in the object graph are implemented by the :py:class:`.OGO` class).
Each of the objects is from one of the :term:`pipeline`'s object types and is assigned with an *object id* that must be unique string 
among all objects for the same object type (i.e. there can be only one object of type *individual* with object id *john*).

Each object is also associated with list of dependency objects. The dependency objects are objects whose targets will be used 
in the creation of the targets of the current object. A target, ``T(t)`` of an object is created by the rule from the snake 
file of the object's object type
that has the target in its output clause (i.e. ``ouput: T(t)``). The input clause of the rule may contain other targets from the 
same object type (:py:func:`.T`), targets in a dependency object (:py:func:`.DT`), or other files.

In addition, each object is associated with project parameters, a dictionary of parameter name to parameter value strings that
provide important information for the creation of the objects targets. 
For example:

.. code-block::

   def run(project,OG):
       params = {"symlink.sample.bai":"<path_bai>","symlink.sample.bam":"<path_bam>"}
       dep = [<list of objects>]
       OG.add(t,i,params,dep)

will create simbolic links ``sample.bam`` and ``sample.bai`` to corresponding paths in the directory for object ``i`` of type ``t``. 

Object are typically created by the *pipelines*'s ``build_object_graph.py`` script with the :py:meth:`~snakeobjects.ObjectGraph.add` method of the :py:class:`.ObjectGraph`.
The order of the dependency objects is preserved and the :py:func:`.DT` and :py:func:`.DP` functions will use the order in the 
bread-first traversal of the object graph. 

During the execution of the *workflow* targets for of the project's object get created and stored
in file in the ``objects`` subdirectory of the *project directory*. 


