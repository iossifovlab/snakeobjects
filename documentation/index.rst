
snakeobjects
============


Overview
--------

**snakeobjects** is a workflow managment framework based on ``snakemake`` that
uses an object-oriented abstraction of workflows. ``snakeobjects`` workflows
are easier to develop, to maintain and to adopt compared to the equivallent
workflows written in ``stanakemake``, but inherit all the powerfull features of
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
using snakemakes syntax where the inputs and outputs specify targets instead of
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

``build_boject_graph.py`` script can use the :py:class:`.Project` to access the
project attributes that point to the relevant meta data and the simple inteface
of :py:class:`.ObjectGraph` to construct the object graph.  

The object-type snakefiles contain
list the targets for the object types and the rules for creating the targets.
The rules are written using the ``snakemake``'s syntax and use a set of extension functions (see
:ref:`snake-extenssions`)  for referring to targets in the current object
(:py:func:`.T`); to targets in dependency objects (:py:func:`.DT`); to
parameters of the current object (:py:func:`.P`); to parameters of the
dependency objects (:py:func:`.DP`); and to global object parameters
(:py:func:`.PP`).  In addition, the pipeline directory can contain scripts,
conda environment definitions or other artefacts used by the pipeline. 

Projects
--------

.. toctree::
    :maxdepth: 2

    installation 
    snakemakeExt
    gettingStarted
    pythonUtils 

.. rstTests (This is a comment)

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
