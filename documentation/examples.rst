.. _examples:

*******************
Additional Examples
*******************

Snakeobjects Paper Example
==========================


This section explains how to load and execute small project described Supplementary of ``snakeobjects`` paper. We assume that you will work or a Linux or Mac. In addition, we assume that you have a conda or miniconda installed (`Conda
Installation
<https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_, to speed up the environment createtion `install Mamba <https://mamba.readthedocs.io/en/latest/installation.html>`_ ).
Everything else needed is inclucded in the
:download:`snakeobjectsPaperExample.tgz <./snakeobjectsPaperExample.tgz>`. When you
download and extract (``tar xzf snakeobjectsPaperExample.tgz``) file, you will
get a directory called ``snakeobjectsPaperExample``. We assume here that this directory is placed in /tmp, but you can place it elsewhere.
Next execute following commands:

.. code-block:: bash
		
	(base) .............$ cd /tmp/snakeobjectsPaperExample/pipeline

	(base) /tmp/snakeobjectsPaperExample/pipeline$ mamba env create

	(base) /tmp/snakeobjectsPaperExample/pipeline$ conda activate snakeobjectsPaperExample

        (snakeobjectsPaperExample) /tmp/snakeobjectsPaperExample/pipeline$ cd ../project
	
	(snakeobjectsPaperExample) /tmp/snakeobjectsPaperExample/project$ sobjects prepare

        (snakeobjectsPaperExample) /tmp/snakeobjectsPaperExample/project$ sobjects run -j -q

In the subdirectory denovos/all you will find allDenovoCalls.txt file.
Opening it in excel with a little formatting can get:

.. image:: _static/paperExample-allDenovoCalls.png
  :alt: allDenovoCalls.txt


Snakemake Tutorial Example
==========================

This section explains how to load and execute the small project reimplementing
in snakeobject
`Snakemake Tutorial example <https://Snakemake.readthedocs.io/en/stable/tutorial/tutorial.html>`_.
We assume that you will work or a Linux or Mac. In addition, we assume that you
have conda or miniconda installed (`Conda
Installation
<https://docs.conda.io/projects/conda/en/latest/user-guide/install>`_).
Everything else needed is inclucded in the
:download:`SnakemakeTutorialExample.tgz <./SnakemakeTutorialExample.tgz>`. When you
download and extract (``tar xzf SnakemakeTutorialExample.tgz``) file, you will
get a directory called ``snakeomakeTutorialExample``. We assume here that this directory is placed in /tmp, but you can place it elsewhere.
You can then follow the same steps as described in the ``Paper Supplementary example`` above.

Tests Demo Examples
===================

There are nine mini projects in the directory ``snakeobjects/tests/demos``.
Each mini project describes some features of snakeobjects.
