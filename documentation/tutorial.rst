********
Tutorial
********

Introduction
============

In this tutorial we will gradually build a complex example that mimics closely a 
realistic scenario using simulated data. We will pretend to be a bioinformatician 
involved in a large project for examining *de novo* mutations in a set of 
children. We were just given a directory  containing the first batch of exome sequence data 
generated from 100 families comprised of mother, father, and a child and we were asked to examine the 
quality of the data and to identify the *de novo* substitutions in the 100 children. A *de novo* 
substitution is a nucleotide at a give position in a child that is not present in the her parents. 

Setup
=====

We assume that the ``snakeobjects`` is already installed following the
installation  instruction in :ref:`getting-started-label`. We also assume that
the active conda environment has pandas and matplotlib installed in addition to
the snakeobjects.  A commands like ``conda install pandas`` and ``conda install
matplotlib`` whould be sufficient.

The contents of the directory we were given can be downloaded from
:download:`snakeobjects_tutorial_input.tgz
<./snakeobjects_tutorial_input.tgz>`. When you download and extract this file
(on linux and Mac that is done by ``tar xzf snakeobjects_tutorial_input.tgz``)
you will get a directory called ``input``. Inside you find a shorted
``README.txt`` file describing the contents. 

Briefly, the directory contains the reference genome (``chrAll.fa``), the known
genes (``genes.txt``), the regions targeted by the exome capture
(``targetRegions.txt``).  For the purposes of the tutorial we will use only the
1MB for the human chromosome 1, that contains 54 transcripts from 36 different
genes, that include 67 protein coding regions targeted by the exome capture.

There are also 764 fastq files in the sub-directory ``fastq`` grouped into 384
pairs of read 1---read 2 files (i.e.  ``fastq/FC0A03F0C/L007/bcL_R1.fastq.gz``
and ``fastq/FC0A03F0C/L007/bcL_R2.fastq.gz``. This a typical structure of the
output for the so-called paired-end sequencing. In addition, there is a file
called ``fastqs.txt`` providing the additional information for each fo the
fastq files (i.e. the individuals for each of the fastq files) and the
``collection.ped`` file describing the 100 trio families using a standard
pedigree file that has one line for each fo the 300 individuals. 

Step 1. Examining the fastq files
=================================

We will build a pipeline for achieving our goals or identifying the *de novo*
variants and the judging the qaulity of the data by small steps. We will start
with a simple pipeline that checks if the fastq files are valid and counts the
number of pairs of reads for each fastq run. 

Step 1.1. Create project and pipline directories
------------------------------------------------

First, let's crete a directory called pipeline where will create add the
pipeline's componenets and a directory called project where will store the
results of the pipeline. 

.. code-block:: bash
    
    $ mkdir pipeline 
    $ mkdir project

Below, we assume that these directories are 'next to' each other, but this is
not required by snakeobjects.  We will aslo assume that the input directory has
been extracted next to the pipeline and project directory.

Step 1.2. Configure the projects
--------------------------------

Next, we will create a file called ``so_project.haml`` in the project directory
with the following content:

.. code-block:: bash
    
    so_pipeline: "[P:projectDir]/../pipeline"

    inputDir: "[P:projectDir]/../input"
    fastqsFile: "[C:inputDir]/fastqs.txt" 
    fastqDir: "[C:inputDir]/fastq" 


This ``so_project.yaml`` configures our first snakeobjects project. The first
line specifies that the pipeline that will operate on this project is contained
in the directory ``../pipeline`` relative to the project directory. The
``[P:projectDir]`` is an example for *interpolation* feature available in the
project configuration files. Alternatively (and possibly preferably), we could
have specified the pipeline directory using a full path. The second line adds a
configuration property called ``inputDir`` that points to the input directory
that was given to us. Here we have also configured relatively to the project
directory which is convenient for the purposes of the tutorial. The
``inputDir`` will not directory be used by the pipeline, but is convenient to
define several of the project properties relative to the ``inputDir``. Indeed,
the next two properties  ``fastqsFile`` and ``fastqDir`` are defined succinctly
based on the ``inputDir`` and point to the table describing the fastq files and
to the directory containing the fastq files. Using ``inputDir`` is optional ---
we could have defined the fastqDir as ``"[P:projectDir/../input/fastq"`` or
with a full path to the directory.

To check if the configuration was successful, we can use the :option:`sobjects describe`
command from within the project directory (the easiest ways to tell ``sobjects`` which project
it should operate on is to 'go into' the project's directory or any of its subdirectories):

.. code-block:: bash

    $ cd project
    $ sobjects describe
    # WORKING ON PROJECT /home/iossifov/work/snakeobjects/tutorial/project
    # WITH PIPELINE /home/iossifov/work/snakeobjects/tutorial/pipeline
    Project parameters:
        so_pipeline: /home/iossifov/work/snakeobjects/tutorial/project/../pipeline
        inputDir: /home/iossifov/work/snakeobjects/tutorial/project/../input
        fastqDir: /home/iossifov/work/snakeobjects/tutorial/project/../input/fastq
        fastqsFile: /home/iossifov/work/snakeobjects/tutorial/project/../input/fastqs.txt

The result should show that sobjects has determined the project and the pipeline directories 
and that the fastqdir and fastqsFile project properties point to the correct locations:

.. code-block:: bash

    $ head /home/iossifov/work/snakeobjects/tutorial/project/../input/fastqs.txt
    flowcell	lane	barcode	individual
    FC0A03F0F	L004	J	SM07279
    FC0A03F0F	L004	K	SM04710
    FC0A03F0F	L004	L	SM63089
    FC0A03F0C	L007	J	SM18469
    FC0B03F00	L001	J	SM18469
    FC0A03F0C	L007	K	SM64466
    FC0B03F00	L001	K	SM64466
    FC0A03F0C	L007	L	SM78901
    FC0B03F00	L001	L	SM78901

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

.. code-block::  

    import pandas as pd
    from pathlib import Path

    def run(proj, OG):
        fastqDir = Path(proj.parameters['fastqDir'])
        fastqs = pd.read_table(proj.parameters["fastqsFile"], sep='\t', header=0)

        for i, r in fastqs.iterrows():
           OG.add('fastq',
                   ".".join([r['flowcell'],r['lane'],r['barcode']]),
                   {
                     'R1': fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R1.fastq.gz",
                     'R2': fastqDir / r['flowcell'] / r['lane'] / f"bc{r['barcode']}_R2.fastq.gz",
                     'sampleId': r['individual'],
                   }
           )

The ``run`` function is given the project (``proj``) for which it will create a
new object graph and object graph instance (``OG``) that add the new
objects into. The function uses proj.parameters to access the necessary
parameters, the ``fastqDir`` pointing to the directory with the fastq files and
the ``fastqsFile`` pointing to the table describing the project fastq pairs of
files, the first few lines of which are shown above. The function, uses the
pandas to read and iterate over all lines of this table and adds (:py:meth:`snakeobjects.ObjectGraph.add`) 
an object of type ``fastq`` and object id equal to
the ``flowcell``, ``lane``, and ``'barcode`` properties concatenated with ``.``. For example, the
object created for the first line of the fastqs.txt file will have an object id
equal to ``FC0A03F0F.L004.J``. Three parameters are also added to each of
objects: ``R1`` and ``R2`` point to the fastq files for the first and for the
second reads defined relative to the project's ``fastqDir`` parameter, and the 
``sampleId`` is assigned the value of the ``individual`` column.

Step 1.4. Prepare the projects 
------------------------------

Next we will create the object graph for our project. We do that by using the :option:`sobjects prepare` command
from within the project directory. We can then follow with the :option:`sobjects describe` to see 
description of the created object graph:

.. code-block:: bash

    $ sobjects prepare
    # WORKING ON PROJECT /home/iossifov/work/snakeobjects/tutorial/project
    # WITH PIPELINE /home/iossifov/work/snakeobjects/tutorial/pipeline-step-1

    $ sobjects describe
    # WORKING ON PROJECT /home/iossifov/work/snakeobjects/tutorial/project
    # WITH PIPELINE /home/iossifov/work/snakeobjects/tutorial/pipeline-step-1
    Project parameters:
        so_pipeline: /home/iossifov/work/snakeobjects/tutorial/project/../pipeline-step-1
        inputDir: /home/iossifov/work/snakeobjects/tutorial/project/../input
        fastqDir: /home/iossifov/work/snakeobjects/tutorial/project/../input/fastq
        fastqsFile: /home/iossifov/work/snakeobjects/tutorial/project/../input/fastqs.txt
    Object types:
         fastq : 384

The above says that we have created an object graph that has 384 objects of type ``fastq``, 
that is exactly what we expected. 

Importantly, the :option:`sobjects prepare` command created a directory ``objects`` in 
the project directory. The ``objects`` contains large number of subdirectories and two files:

.. code-block:: bash

    $ find objects | head
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

.. code-block::
   
    add_targets() 

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

    $ sobjects run -j -q
    # WORKING ON PROJECT /mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/project
    # WITH PIPELINE /mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/project
    export SO_PIPELINE=/mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/project/objects/.snakeobjects/main.snakefile -d /mnt/wigtop1/home/iossifov/work/snakeobjects/tutorial/project/objects -j -q
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
``so_fastq_ojb`` rule is used to build an object of type ``fastq`` and since we
have 384 such objects in our graph this rule is used for 384 jobs.  

The :option:`sobjects run` command added a ``objects/.snakemake`` directory
where ``snakemake`` stored internal information related to the execution. This
may come handy for figuring out error in complex situation but we will no 
cover the ``snakemake`` privates in this tutoria and refer you to the ``snakemake``'s 
documentation for more information. 
Most importantly, :option:`sobjects run` created a ``obj.flag`` file in directories 
for each object: 

.. code-block:: bash

    $ find objects/fastq  | head
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
soluction-step-1.5 directory.

Step 1.6. Add a usefull target 
------------------------------

Next, we will add an explicit target that does something useful: we will create a target
called ``pairNumber.txt`` for objects of type ``fastq`` that will store the number of paired reads
in the fastq object. This is achieved by replacing the auto-generated one line in the pipeline's 
``fastq.snakemake`` with the following:

.. code-block:: 
    :linenos:

    add_targets("pairNumber.txt")

    rule countReads:
        input: P('R1'), P('R2')
        output: T("pairNumber.txt")
        run:
            import gzip

            nPairs= 0
            buff = []
            with gzip.open(input[0]) as R1F, \
                 gzip.open(input[1]) as R2F:
                for l1,l2 in zip(R1F,R2F):
                    buff.append((l1,l2))
                    if len(buff) == 4:
                        nPairs += 1
                        buff = []

            with open(output[0],"w") as OF:
                OF.write(f'{wildcards.oid}\t{nPairs}\n')


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

We use the ``run`` clause of the ``countReads`` to implement the generation of
the outputs from inputs by a python snipped. There are alternative clauses that
can be used instead of ``run`` that allow different ways to implement the
output generation (i.e. whit ``shell`` clause we can use shell commands,
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

Step 1.7. Crate a test project
------------------------------

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

    $ mkdir projectTest
    $ cd projectTest
    $ cp ../project/so_project.yaml
    $ head -10 ../input/fastqs.txt  > fastqs-small.txt

and replace line 5 in ``projectTest/so_project.yaml`` file
that configures the ``fastqsFile`` project parameter to
point to the ``fastqs-small.txt`` file containing a description of 9 fastq runs
instead of the complete ``fastqs.txt`` file with 384 fastq runs:

.. code-block::
    :linenos:
    :emphasize-lines: 5

    so_pipeline: "[P:projectDir]/../pipeline"

    inputDir: "[P:projectDir]/../input"
    fastqDir: "[C:inputDir]/fastq"
    fastqsFile: "[P:projectDir]/fastqs-small.txt"

With the ``projectTest`` configured, we can then *prepare* and *run* the projects:

.. code-block:: bash

    $ sobjects run -j -q
    # WORKING ON PROJECT /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/projectTest
    # WITH PIPELINE /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/pipeline

    $ sobjects run -j -q
    # WORKING ON PROJECT /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/projectTest
    # WITH PIPELINE /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/pipeline
    UPDATING ENVIRONMENT:
    export SO_PROJECT=/Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/projectTest
    export SO_PIPELINE=/Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/pipeline
    export PATH=$SO_PIPELINE:$PATH
    RUNNING: snakemake -s /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/projectTest/objects/.snakeobjects/main.snakefile -d /Users/iiossifov/work/snakeobjects/tutorial/solution-step-1.8/projectTest/objects -j -q
    Job counts:
        count	jobs
        9	countReads
        1	so_all_targets
        9	so_fastq_obj
        19
    ...

The run finishes almost instantaneously and as a results we can find the pairNumber.txt files 
for each of the 9 fastq objects created for the projectTest:

.. code-block:: bash

    $ cat objects/fastq/*/pairNumber.txt 
    FC0A03F0C.L007.J	2038
    FC0A03F0C.L007.K	1994
    FC0A03F0C.L007.L	2086
    FC0A03F0F.L004.J	2238
    FC0A03F0F.L004.K	2212
    FC0A03F0F.L004.L	2278
    FC0B03F00.L001.J	214
    FC0B03F00.L001.K	245
    FC0B03F00.L001.L	223

We can now spot check to see if the reported number of reads is correct. For example: 

.. code-block:: bash

    $ cat ../input/fastq/FC0A03F0C/L007/bcJ_R1.fastq.gz  | gunzip -c | wc -l
    8152

Read 1 file for the fastq run ``FC0A03F0C.L007.J`` contains 8152 lines which is equal to 4 times 
the number of pairs (2038) reported in the corresponding pairNumber.txt file. This is exactly what is 
expected: as described above, sequencing reads are represented in 4 lines in the fastq files. 

Step 1.8. Re-run the project 
----------------------------

Note. If we had a cluster profile configured we can 
used it!!!

Step 1.9. Add a summary object 
------------------------------

draw a graph??

