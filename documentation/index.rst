
snakeobjects
============

.. image:: https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod
    :target: https://gitpod.io/#https://github.com/snakemake/snakemake

.. image:: https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/github/workflow/status/snakemake/snakemake/Publish%20to%20Docker%20Hub?color=blue&label=docker%20container&branch=master
    :target: https://hub.docker.com/r/snakemake/snakemake

.. image:: https://github.com/snakemake/snakemake/workflows/CI/badge.svg?branch=master&label=tests
    :target: https://github.com/snakemake/snakemake/actions?query=branch%3Amaster+workflow%3ACI

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: https://stackoverflow.com/questions/tagged/snakemake

.. image:: https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow
    :target: https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester

.. image:: https://img.shields.io/discord/753690260830945390?label=discord%20chat   
    :alt: Discord
    :target: https://discord.gg/NUdMtmr

.. image:: https://img.shields.io/github/stars/snakemake/snakemake?style=social
    :alt: GitHub stars
    :target: https://github.com/snakemake/snakemake/stargazers

.. .. raw:: html
          <span class="__dimensions_badge_embed__" data-doi="https://doi.org/10.1093/bioinformatics/bts480" data-legend="always" data-style="large_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

Overview
--------

**snakeobjects** is a workflow management framework based on ``snakemake`` that
uses an object-oriented abstraction of workflows. ``snakeobjects`` workflows
are easier to develop, to maintain and to adopt compared to the equivallent
workflows written in ``stanakemake``, but inherit all the powerful features of
``snakemake``. These include the portability, efficient resource usage, the
large expressive power due to the tight python integration, and the large
community of the ``snakemake`` users.

.. toctree::
    :maxdepth: 1

    gettingStarted
    tutorial
    examples
    documentation
    specialTopics
    terms
    sobjects
    snakemakeExt
    pythonUtils

.. rstTests (This is a comment)

Documentation
-------------

https://www.iossifovlab.com/snakeobjects

Source
------

https://github.com/iossifovlab/snakeobjects

Setting up the development environment
--------------------------------------

.. code-block:: bash

    git clone git@github.com:iossifovlab/snakeobjects.git
    cd snakeobjects
    conda env create environment.yml
    conda activate snakeobjectsDev
    pip install -e .

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
