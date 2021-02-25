********
Tutorial
********

Introduction
============

In this tutorial we will gradually build a complex example that mimics closely
a realistic scenario using simulated data. We will pretend to be a
bioinformatician involved in a large project for examining *de novo* mutations
in a set of children. We were just given a directory  containing the first
batch of exome sequence data generated from 100 families comprised of mother,
father, and a child and we were asked to examine the quality of the data and to
identify the *de novo* substitutions in the 100 children. A *de novo*
substitution is a nucleotide at a given position in a child that is not present
in her parents. 

Setup
=====

``snakeobjects`` and this tutorial are tested and work well on Linux; they
don't work on Windows. So we assume that you will work on a Linux on Mac. In
addition, we assume that you have a conda or miniconda installed (`Conda
Installation
<https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_).
Everything else needed to follow this tutorial is inclucded in the
:download:`snakeobjectsTutorial.tgz <./snakeobjectsTutorial.tgz>`. When you
download and extract (``tar xzf snakeobjectsTutorial.tgz``) the file, you will
get a directory called ``snakeobjectsTutorial``. In the examples that follow we
have assumed that the ``snakeobjectsTutorial`` directory is created in the
``/tmp`` directory, but this is not essential. Everything will work just fine
if you have the ``snakeobjectsTutorial`` in your home directory or anywhere
else.

The ``snakeobjectsTutorial`` folder contains one file named ``environment.yml``
and two subdirectories: ``input``, and ``solutions``.

The ``envirnoment.yml`` file defines the conda environment we will use
troughout the tutorial:

 
.. literalinclude:: snakeobjectsTutorial/environment.yml

To create the ``snakeobjectsTutorial`` evironment and to activate it we can use
the follwoing commands:

.. code-block:: bash

    (base) .........................$ cd /tmp/snakeobjectsTutorial 
    (base) /tmp/snakeobjectsTutorial$ conda env create 
    (base) /tmp/snakeobjectsTutorial$ conda activate snake 
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ 

The command line prompt shows the activate conda environment and the current
working directory. We will show this long prompt every where below to remind
you that you must work under the ``snakeobjectsTutorial`` environment and
because a lot of the commands depend on the current working directory.  

The ``input`` directory naturally contains the input data for our project.
Inside it you will find a short ``README.txt`` file describing the contents.
Briefly, the directory contains the reference genome (``chrAll.fa``), the known
genes (``genes.txt``), and the regions targeted by the exome capture
(``targetRegions.txt``).  For the purposes of the tutorial we will use only the
first 1MB of the human chromosome 1 that contains 54 transcripts from 36
different genes. The genes include 67 protein coding regions/exons that are
targeted by the hypothetical exome capture. 
There are also 768 fastq files in the sub-directory ``input/fastq`` grouped
into 384 pairs of read 1---read 2 files (i.e.
``fastq/FC0A03F0C/L007/bcL_R1.fastq.gz`` and
``fastq/FC0A03F0C/L007/bcL_R2.fastq.gz``. This a typical structure of the
output for the so-called paired-end sequencing. In addition, there is a file
called ``fastqs.txt`` providing the additional information for each of the
fastq files (i.e. the individuals for each of the fastq files) and the
``collection.ped`` file describing the 100 trio families using a standard
pedigree file that has one line for each of the 300 individuals. 

Finally, the ``solutions`` subdirectory contains the state of the pipeline and
the projects at various stages (steps) of the tutorial. For example,
``solutions/step-2.4`` contains the state after we have completed Step 2.4. For
the impatient, the ``solutions/final`` contains the complete pipeline built by
the end of the tutorial. 

Step 1. Examining the fastq files
=================================

We will build a pipeline for achieving our goals or identifying the *de novo*
variants and the judging the quality of the data by small steps. We will start
with a simple pipeline that counts the number of pairs of reads for each fastq run. 

Step 1.1. Create project and pipline directories
------------------------------------------------

First, let's create a directory called ``pipeline`` where will create add the
pipeline's components and a directory called ``project`` where will store the
results of the pipeline. We will assume that you create these directories in 
the ``snakeobjectsTutorial`` directory:  

.. code-block:: bash
    
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ mkdir pipeline 
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ mkdir project 

This location of the directories is not required. Provided you are willing
to make few minor changes in the projects' configuration you can create
the two directories anywhere. 

Step 1.2. Configure the projects
--------------------------------

Next, we will create a file called ``so_project.yaml`` in the ``project`` directory
with the following content:


.. literalinclude:: snakeobjectsTutorial/solutions/step-1.5/project/so_project.yaml 
    

This ``so_project.yaml`` configures our first snakeobjects project. The first
line specifies that the pipeline that will operate on this project is contained
in the directory ``../pipeline`` relative to the project directory. The
``[P:projectDir]`` is an example for *interpolation* feature available in the
project configuration files. Alternatively (and possibly preferably), we could
have specified the pipeline directory using a full path (``/tmp/snakeobjectsTutorial/pipeline``). The second line adds a
configuration property called ``inputDir`` that points to the input directory
that was given to us. Here we have also configured relatively to the project
directory which is convenient for the purposes of the tutorial. The
``inputDir`` will not directory be used by the pipeline, but is convenient to
define several of the project properties relative to the ``inputDir``. Indeed,
the next two properties  ``fastqsFile`` and ``fastqDir`` are defined succinctly
based on the ``inputDir`` and point to the table describing the fastq files and
to the directory containing the fastq files. Using ``inputDir`` is optional ---
we could have defined the fastqDir as ``"[P:projectDir/../input/fastq"`` or
with a full path to the directory, ``/tmp/snakeobjectsTutorial/input/fastq``.

To check if the configuration was successful, we can use the :option:`sobjects describe`
command from within the project directory (the easiest way to tell ``sobjects`` which project
it should operate on is to 'go into' the project's directory or any of its subdirectories):

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ $ cd project
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects describe
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    Project parameters:
        so_pipeline: /tmp/snakeobjectsTutorial/project/../pipeline
        inputDir: /tmp/snakeobjectsTutorial/project/../input
        fastqDir: /tmp/snakeobjectsTutorial/project/../input/fastq
        fastqsFile: /tmp/snakeobjectsTutorial/project/../input/fastqs.txt
    Object types:


The result should show that sobjects has determined the project and the pipeline directories 
and that the fastqdir and fastqsFile project properties point to the correct locations:

.. code-block:: bash

    ....$ head /tmp/snakeobjectsTutorial/project/../input/fastqs.txt
    flowcell	lane	barcode	individual
    FC0A03F09	L006	D	SM90370
    FC0B03F03	L001	D	SM90370
    FC0A03F09	L006	E	SM03231
    FC0B03F03	L001	E	SM03231
    FC0A03F09	L006	F	SM79279
    FC0B03F03	L001	F	SM79279
    FC0A03F0C	L007	D	SM14701
    FC0B03F00	L008	D	SM14701
    FC0C03F00	L004	D	SM14701


Step 1.3. Create the build_object_graph.py 
------------------------------------------

Now that have configured our first project, we will turn our attention to the 
pipeline. So far the pipeline directory is empty. The first thing to do when 
starting a pipeline is to create the ``build_object_graph.py`` script. In the 
**Step 1** we will create a very simple graph that contains one object for each 
fastq pairs of files. The fastq pairs are listed in the fastqs.txt file in the input
directory and we have already ensured that our project has a prarameters 
(``fastqsFile``) that points on the fastqs.txt. The contents of our first 
``build_object_graph.py`` are shown below. You should create a file named 
``build_object_graph.py`` in the pipeline directory and copy the shown contents
in the file.

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.5/pipeline/build_object_graph.py

The ``run`` function is given the project (``proj``) for which it will create a
new object graph and object graph instance (``OG``) that it should add the new
objects into. In the first two lines, the function accesses the necessary
parameters through ``proj.parameters``: the ``fastqDir`` pointing to the
directory with the fastq files and the ``fastqsFile`` pointing to the table
describing the project fastq pairs of files (or fastq runs).  The function,
then uses pandas to read and iterate over all lines of this table and for each
line adds (:py:meth:`snakeobjects.ObjectGraph.add`) an object of type ``fastq``
and object id equal to the ``flowcell``, ``lane``, and ``'barcode`` properties
concatenated with ``.``. For example, the object created for the first line of
the fastqs.txt file will have an object id equal to ``FC0A03F0F.L004.J``. Three
parameters are also added to each of objects: ``R1`` and ``R2`` point to the
fastq files for the first and for the second reads defined relative to the
project's ``fastqDir`` parameter, and the ``sampleId`` is assigned the value of
the ``individual`` column.

Step 1.4. Prepare the projects 
------------------------------

Next we will create the object graph for our project. We do that by using the :option:`sobjects prepare` command
from within the project directory. We can then follow with the :option:`sobjects describe` to see 
description of the created object graph:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects prepare
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects describe
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    Project parameters:
        so_pipeline: /tmp/snakeobjectsTutorial/project/../pipeline
        inputDir: /tmp/snakeobjectsTutorial/project/../input
        fastqDir: /tmp/snakeobjectsTutorial/project/../input/fastq
        fastqsFile: /tmp/snakeobjectsTutorial/project/../input/fastqs.txt
    Object types:
         fastq : 384

The above says that we have created an object graph that has 384 objects of type ``fastq``, 
that is exactly what we expected. 

Importantly, the :option:`sobjects prepare` command created a directory ``objects`` in 
the project directory. The ``objects`` contains large number of subdirectories and two files:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ find objects | head
    objects
    objects/.snakeobjects
    objects/.snakeobjects/OG.json
    objects/.snakeobjects/main.snakefile
    objects/fastq
    objects/fastq/FC0A03F09.L006.G
    objects/fastq/FC0A03F09.L006.G/log
    objects/fastq/FC0A03F0F.L007.E
    objects/fastq/FC0A03F0F.L007.E/log
    objects/fastq/FC0A03F06.L008.I

The ``objects/.snakeobjects/OG.json`` file stores the object graph that was just created
and the ``objects/.snakeobjects/main.snakefile`` is the projects specific snakefile that 
will be provided to the snakemake upon execution of the pipeline. In addition, there are
directories for each object from the objects graph where the objects' targets will be 
stored: the ``objects/fastq/FC0A03F09.L006.G`` directory will contain the targets
for the object of type ``fastq`` and object id ``FC0A03F09.L006.G``. Each of the 
object directories has also a log subdirectory (i.e. ``objects/fastq/FC0A03F09.L006.G/log``) 
where log files associated with the object will be stored (more about log files later). 

In addition, the :option:`sobjects prepare` created one file, ``fastq.snakefile``, in 
pipeline directory. This is not a typical behaviour: as a rule **sobjects** only updates
the project directory (and to be more specific, only its ``objects`` subdirectory),
but when we start a new pipeline it's handy to have placeholders for the object 
type snakemake files be created for us. The content if the new ``fastq.snakefile`` is
very simple: 

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.5/pipeline/fastq.snakefile

This simple one line accomplishes nothing but reminding us that next step would be to
declare the targets to be created for objects of  type ``fastq``. There is one special
target automatically added to every object file, and it is called ``T("obj.flag")``. It is
created after all the other targets for the object are successfully created. Without 
explicitly adding object type targets (as in the current state of the ``fastq.snakefile``), 
the ``T("obj.flag")`` is the only target. We will keep this situation for now, and add 
useful targets shortly.

Step 1.5. Execute the dummy project 
-----------------------------------

After we have *prepared* the project, it time to execute it. This is done 
by :option:`sobjects run` command. This commands requires at least one argument to 
that we use to control how to execute the pipeline:  ``-j`` means that we can use all 
the available local cores to execute the pipeline; ``-j 1`` means that we can use only 
one of the local cores; ``-j 2`` means that we can use 2 local cores, etc. Alternatively, 
we can add ``--profile <cluster profile>`` option which indicate that we will execute 
the pipeline by submitting jobs to the computation cluster defined by the ``<cluster profile>``.
Executing pipelines on cluster will be covered later. For now we would use the simplest 
``-j`` form to execute the *dummy* pipeline we have developed so far. We can specify more options
to the :option:`sobjects run` command that control the execution of the pipeline and here we would
use the ``-q`` command that instructs the underlying ``snakemake`` to be 'quiet' and produce only 
minimal output:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects run -j -q
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/project
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/project/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/project/objects -j -q
    Job counts:
        count	jobs
        1	so_all_targets
        384	so_fastq_obj
        385

As usual, the first two lines show the project and the pipeline that
``sobjects`` operates with. The next five lines provide information about how
``sobjects`` executes ``snakemake``, including the environment variables and
the command line parameters passed to ``snakemake``. Finally, we see the number
of successfully executed jobs summarized  by the ``snakemake`` rules used for
each job.  Both rules shown are automatically generated by ``snakeobjects`` as
indicated by the ``so_`` prefix. The ``so_all_targets`` is the default (first)
rule ``snakemake`` sees and is the one that specifies all targets that need to
be created for the project and is naturally executed only 1 time. The
``so_fastq_obj`` rule is used to build an object of type ``fastq`` and since we
have 384 such objects in our graph this rule is used for 384 jobs.  

The :option:`sobjects run` command added a ``objects/.snakemake`` directory
where ``snakemake`` stored internal information related to the execution. This
may come handy for figuring out errors in complex situation but we will not 
cover the ``snakemake`` privates in this tutoria and refer you to the ``snakemake``'s 
documentation for more information. 
Most importantly, :option:`sobjects run` created a ``obj.flag`` file in directories 
for each object: 

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ find objects/fastq  | head
    objects/fastq
    objects/fastq/FC0A03F09.L006.G
    objects/fastq/FC0A03F09.L006.G/obj.flag
    objects/fastq/FC0A03F09.L006.G/log
    objects/fastq/FC0A03F0F.L007.E
    objects/fastq/FC0A03F0F.L007.E/obj.flag
    objects/fastq/FC0A03F0F.L007.E/log
    objects/fastq/FC0A03F06.L008.I
    objects/fastq/FC0A03F06.L008.I/obj.flag
    objects/fastq/FC0A03F06.L008.I/log


These files correspond the ``T("obj.flag")`` target and are created only after all
of the other targets for the object are created. As we have not yet added any other 
targets, only the ``obj.flags`` are create for the objects. 

The pipline and the project configuration we have developped so far are included in the 
``solutions/step-1.5 directory``.

Step 1.6. Add a usefull target 
------------------------------

Next, we will add an explicit target that does something useful: we will create a target
called ``pairNumber.txt`` for objects of type ``fastq`` that will store the number of paired reads
in the fastq object. This is achieved by replacing the auto-generated one line in the pipeline's 
``fastq.snakemake`` with the following:

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.8/pipeline/fastq.snakefile
    :linenos:

The first line declares that objects of type ``fastq`` (the object type is
implied by the name of the snakemake file ``fastq.snakemake``) will have a
target named ``pairNumber.txt``.

Next is a ``snakemake`` rule we have called ``countReads`` that describes how
such target is to be created: the output clause on line 5 shows that this rule
will generate the ``T("pairNumber.txt")`` target.  As input (line 4), the rule
will use the two files pointed by the fastq object's parameters called ``R1``
and ``R2``.  To obtain the value of ``R1`` and ``R2`` parameters for the object
the rule operates on the rule uses the extension function :py:func:`.P`. As a
reminder, these parameters were set so that their values are full file paths to
the read 1 and read 2 files object within the ``input/fastq`` directory for the
fastq object.

We use the ``run`` clause of the ``countReads`` rule to implement the generation of
the outputs from inputs by a python snipped. There are alternative clauses that
can be used instead of ``run`` that allow different ways to implement the
output generation (i.e. with ``shell`` clause we can use shell commands,
like ``cat``, ``echo`` to generate the outputs). Later in the tutorial we will
demonstrate how to use some of these alternatives. 

The python snipped uses the ``input`` and ``output`` objects provided by
``snakemake`` to access the rule's input, outputs. The snipped also uses the
``snakemake``'s object called ``wildcards`` object to access the object id
(``oid``) for the object the rule operates on. (For the readers that are
familiar with ``snakemake``, the ``oid`` can be explained by the fact that the
``output: T("pairNumber.txt")`` is equivalent to ``output:
fastq/{oid}/pairNumber.txt``.)

The actual implementation is fairly trivial assuming one knows a bit about the
way pair-end sequencing results are represented in the fastq files. Briefly, in
pair-end sequencing a small (i.e. 200-500 base-pairs) linear DNA fragments from the
genome are sequenced and for each fragment the sequencing machine first reads
several (i.e. 100) nucleotides from one side (read 1) of the fragment and next
several nucleotides from the other end (read 2) of the fragment.  The read 1s
from all fragments are stored in read 1 fastq file and the read 2s are
stored in the read 2 fastq file.  Importantly, the order of the fragments in
the read 1 and read 2 fastq files is the same, so that the first read in read 1
fastq files is from the same DNA fragment as the first read in the read 2 fastq
file. Reads in a fastq file are represented by 4 consecutive lines (fragment id
line, sequence line, separator line, base-quality line). 

Step 1.7. Create a test project
-------------------------------

We just improved our pipleine to do something usefull. In a large project usefull 
tasks usually take a lot of computational resoruces and a lot or time. To be 
able to easily test if the updated pipeline works well it is a good idea to 
have a small test project configured to operate on a small subset of the input data.
Here we will do just that. Even though the complete tutorial input data is small 
enough that it can be processed on a single processor for less than a minute,
this is a useful demonstration of how easy it is to maintain and operate on 
multiple projects with the same ``snakeobjects`` pipeline. 

We can create the new projectTest with the following simple commands:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ cd .. 
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ mkdir projectTest
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ cd projectTest
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ cp ../project/so_project.yaml
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial$ head -16 ../input/fastqs.txt  > fastqs-small.txt

and replace line 5 in ``projectTest/so_project.yaml`` file
that configures the ``fastqsFile`` project parameter to
point to the ``fastqs-small.txt`` file containing a description of 15 fastq runs
instead of the complete ``fastqs.txt`` file with 384 fastq runs:

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.8/projectTest/so_project.yaml
    :linenos:
    :emphasize-lines: 5

With the ``projectTest`` configured, we can then *prepare* and *run* the projects:

.. code-block:: bash

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/projectTest$ sobjects prepare
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/projectTest
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    (snakeobjectsDev) /tmp/snakeobjectsTutorial/projectTest$ sobjects run -j -q
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/projectTest
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/projectTest
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/projectTest/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/projectTest/objects -j -q
    Job counts:
        count	jobs
        15	countReads
        1	so_all_targets
        15	so_fastq_obj
        31
    ...

The run finishes almost instantaneously and as a results we can find the pairNumber.txt files 
for each of the 9 fastq objects created for the projectTest:

.. code-block:: bash

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/projectTest$ cat objects/fastq/*/pairNumber.txt
    FC0A03F09.L006.D	942
    FC0A03F09.L006.E	1037
    FC0A03F09.L006.F	1048
    FC0A03F0C.L007.D	1179
    FC0A03F0C.L007.E	1133
    FC0A03F0C.L007.F	1206
    FC0B03F00.L008.D	511
    FC0B03F00.L008.E	483
    FC0B03F00.L008.F	502
    FC0B03F03.L001.D	1205
    FC0B03F03.L001.E	1290
    FC0B03F03.L001.F	1251
    FC0C03F00.L004.D	684
    FC0C03F00.L004.E	614
    FC0C03F00.L004.F	647


We can spot check to see if the reported number of reads is correct. For example: 

.. code-block:: bash

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/projectTest$ cat ../input/fastq/FC0A03F09/L006/bcD_R1.fastq.gz | gunzip -c | wc -l
    3768

Read 1 file for the fastq run ``FC0A03F09.L006.D`` contains 3,768 lines which
is equal to 4 times the number of pairs (942) reported in the corresponding
pairNumber.txt file. This is exactly what is expected: as described above,
sequencing reads are represented in 4 lines in the fastq files. 

Step 1.8. Re-run the project 
----------------------------

Now that have verified that updated pipeline works, it is time to count the
pair numbers for the complete project:

.. code-block:: bash

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/projectTest$ cd ../project
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects run -j -q
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/project
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/project/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/project/objects -j -q
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ cat objects/fastq/*/pairNumber.txt | head
    cat: 'objects/fastq/*/pairNumber.txt': No such file or directory

The results seem strange. ``snakemake`` doesn't seem to run any jobs and the
``pairNumber.txt`` targets are not created. One way to figure out what's going
is to run ``snakemake`` in verbose mode, by removing the ``-q`` flag:
 
.. code-block:: bash
    :emphasize-lines: 10 

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/project$ sobjects run -j 
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/project
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/project/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/project/objects -j
    Building DAG of jobs...
    Nothing to be done.
    Complete log: /tmp/snakeobjectsTutorial/project/objects/.snakemake/log/2021-02-25T114732.563671.snakemake.log

``snakemake`` says that there is nothing to do. This is a peculiar behaviour of
``snakemake`` that manifests anytime we add a new targets and rules to a
pipeline and we need to create the new targets in a projects that was
successfully executed prior to the addition. The way to overcome this is to
force ``snakemake`` to created all targets built by the new rule:

.. code-block:: bash

    (snakeobjectsDev) /tmp/snakeobjectsTutorial/project$ sobjects run -j -R countReads 
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/project
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/project/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/project/objects -j -R countReads
    Building DAG of jobs...
    Using shell: /bin/bash
    Provided cores: 192
    Rules claiming more threads will be scaled down.
    Job counts:
        count   jobs
        384 countReads
        1   so_all_targets
        384 so_fastq_obj
        769
    Select jobs to execute...
    ...
    
This seems to have done the job and now we do indeed have the ``pairCount.txt`` targets:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ cat objects/fastq/*/pairNumber.txt | head
    FC0A03F00.L001.D	2265
    FC0A03F00.L001.E	2181
    FC0A03F00.L001.F	2221
    FC0A03F00.L001.G	2265
    FC0A03F00.L001.H	2170
    FC0A03F00.L001.I	2162
    FC0A03F00.L002.A	1371
    FC0A03F00.L002.B	1276
    FC0A03F00.L002.C	1365
    FC0A03F00.L002.D	1760
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ cat objects/fastq/*/pairNumber.txt | wc -l 
    384


Step 1.9. Add a summary object 
------------------------------

It will be convenient to combine all the ``pairNumber.txt`` files into 
one file that can easily be open in tools like Excel to get a global understanding
of the pair counts across all the fastq runs. Such aggregation need is fairly typical
in workflows and is easy to implement in ``snakeobjects``. We will add one 
object in the object graph that will be responsible for aggregating information 
for all the ``fastq`` objects. We typically use ``xxxSummayry`` as object type
and ``o`` for the object id of such singleton objects. To add our ``fastqSummary/o`` 
objects we add one line to the ``build_object_graph.py`` pipeline script:

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.9/pipeline/build_object_graph.py
    :emphasize-lines: 17 
    
Importantly, we make the ``fastqSummary/o`` object to be dependent by all ``fastq`` objects
already added to the graph, using the ``deps=OG['fastq']`` parameter. 
As descried in the :py:meth:`snakeobjects.ObjectGraph.__getitem__`  method,
indexing an objects graph with an object type returns the list of the objects
of that type that are in the graph. 

We will add two targets to the new ``fastqSummary`` object type by 
creating the ``fastqSummary.snakefile`` file in the pipeline directory with the following 
content:

.. literalinclude:: snakeobjectsTutorial/solutions/step-1.9/pipeline/fastqSummary.snakefile

The first target, called  ``allPairNumers.txt``, is created by the ``gatherPairNumbers``
rule, as indicated by the ``output: T("allPairNumbers.txt")`` clause. This rule 
for the first time in the tutorial shows the use of the most interesting of the 
``snakeojbects``'s functions, :py:func:`.DT`. The :py:func:`.DT` function here specifies
that the input of the rule will be the ``pairNumber.txt`` targets of the objects the
current object depends on. The current objects is of type ``fastqSummary`` (this is clear 
by name of the snakefile) and the object graphs generated by our pipeline here is only 
one object of this type. We we made this object to be dependent on the all the ``fastq`` 
objects. Thus, the input for this rule will be all the ``pairNubmer.txt`` targets of the
``fastq`` objects. The implementation of the rule uses the ``cat`` and ``sort`` unix tools
through the ``shell`` clause. The ``{input}`` in the implementation is replaced by the list
of all ``pairNumber.txt`` files related to the ``pairNumber.txt`` targets for the ``fastq`` objects
and the ``{output}`` is replaced by the single ``allPariNumbers.txt`` output file.

The second target, ``pairNumber.png`` is created by the second rule. This target uses the 
first target as an input and so will be created only after the first one is successfully built.
The implementation is run clause (a python snipped) that uses a ``matplotlib`` to 
draw all a simple graph showing all the pair numbers across the fastq objects. 

Next, we will re-run out two projects. We will redo the ``projectTest`` from scratch to 
make sure that the complete pipeline functions properly:

.. code-block:: bash
    :emphasize-lines: 2 

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ cd ../projectTest
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/projectTest$ rm -r objects
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/projectTest$ sobjects prepare
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/projectTest
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/projectTest$ sobjects run -j -q
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/projectTest
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/projectTest
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/projectTest/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/projectTest/objects -j -q
    Job counts:
        count	jobs
        15	countReads
        1	gatherPairNumbers
        1	pairNumberFigure
        1	so_all_targets
        1	so_fastqSummary_obj
        15	so_fastq_obj
        34
    ....

In the highlighted command we removed the ``objects`` subdirectory and with that we lost all previously 
built targets for the ``projectTest``. So the ``sobject run`` command than needed recreate all
the targets.  For the large ``project`` we will only create the new targets:

.. code-block:: bash

    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/projectTest$ cd ../project
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects prepare
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    (snakeobjectsTutorial) /tmp/snakeobjectsTutorial/project$ sobjects run -j -q
    # WORKING ON PROJECT /tmp/snakeobjectsTutorial/project
    # WITH PIPELINE /tmp/snakeobjectsTutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/tmp/snakeobjectsTutorial/project
    export SO_PIPELINE=/tmp/snakeobjectsTutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /tmp/snakeobjectsTutorial/project/objects/.snakeobjects/main.snakefile -d /tmp/snakeobjectsTutorial/project/objects -j -q
    Job counts:
        count	jobs
        1	gatherPairNumbers
        1	pairNumberFigure
        1	so_all_targets
        1	so_fastqSummary_obj
        4
        ....

Here we did not removed the ``objects`` subdirectory and the targets that were built the last time
we executed the pipeline for the large ``project`` were preserved. The ``sobjects prepare`` command created one new object 
(the ``fastqSummary/o``) and the ``sobjects run`` created its targets. 
Both the ``projectTest`` and ``project`` projects now have an aggregated ``allPairNumbers.txt`` file 
(``.../objects/fastqSummary/o/allPairNumbers.txt``) 
and a ``pairNumber.png`` figure.
(``.../objects/fastqSummary/o/pairNumber.png``). The ``pairNumber.png`` for ``projctTest`` should look like:

.. image:: _static/projectTest-pairNumber.png
  :width: 400
  :alt: pairNumber.png for the projectTest 

and the ``pairNumber.png`` for the large ``project`` should look like:

.. image:: _static/project-pairNumber.png
  :width: 4000 
  :alt: pairNumber.png for the project

The large figure is probably not ideal and its layout can be improved. But the birds eye view is informative
and its resolution is high enough that one can read zoom in enough to see the details.


Step 2. Alignment and target coverage by sample 
===============================================

Step 2.1. Refernce genome indexing
----------------------------------

Step 2.2. Alignment of fastq 
----------------------------

Step 2.3. Merging alignments by individual 
------------------------------------------

Step 2.4. Target coverage by sample and globally 
------------------------------------------------


Step 3. Calling de novo variants
================================

