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

To check if the configuration was successfull, we can use the :option:`sobjects describe`
command from within the project directory:

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
from within the project directory. We can then flow with the :option:`sobjects describe` to see 
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
pipeline directory. This is not a typical behaviour: as a rule sobjects only updates
the project directory (and to be more specific, only its ``objects`` subdirectory),
but when we start a new pipeline it's handy to have placeholders for the object 
type snakemake files place be created for us. The content if the new ``fastq.snakefile`` is
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

Step 1.6. Add a usefull target 
------------------------------

Count number read pairs.

Step 1.7. Crate a test project
------------------------------

Step 1.8. Re-run the project 
----------------------------

Note. If we had a cluster profile configured we can 
used it!!!

Step 1.9. Add a summary object 
------------------------------



