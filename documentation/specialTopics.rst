**************
Special topics
**************

Wokring with clusters
=====================

``Snakeobjects`` provides transparent access to different cluster architectures such as ``slurm``, ``sge``, and others. This is accomplished by placing cluster profile path into ``so_project.yaml`` file with the name ``default_snakemake_args`` as in this example:

.. code-block::

   default_snakemake_args: --profile <path to profile folder>

Profile folder should contain file ``config.yaml`` with directives specific to cluster architecture. For more information on profiles see https://github.com/snakemake-profiles/doc.
Access to profile is accomplished by command ``submit`` ( :option:`sobjects submit` ).



TODO

Multi-part targets
==================

Sometimes target computing time exceeds several hours or even days. This negatively affects overall project management: instabily of file system, failure of computational node, or many other causes may delay processing downstream targets. With multiple processing units it may be valuable to split a target into subtargets. In bioinformatics it is convenient to restrict bam file processing by extracting data for individual chromosomes and then accumulating partial results in a separate target. In ``snakeobjects`` we have special functions that make the subdivision of a target into smaller ones with the following aggregation of these parts into original target easy to accomplish. Below is an example.

.. code-block:: python

    ids = list_of_part_ids

    rule beg:
	  output:
        T('beg-{c}.txt')
      shell:
        "initialize.sh {wildcards.oid} {wildcards.c} > {output}"

    rule part:
      input:
        T('beg-{c}.txt')
      output:
        T('part-{c}.txt')
      shell:
        "process_part.sh {wildcards.oid} {wildcards.c} > {output}; "

    rule mergedI:
      input:
        expand(TE('part-{c}.txt'),c=ids)
      output:
        T('merged.txt')
      shell: 
        "merge_parts.sh {input} > {output}"

Here initialize.sh, process_part.sh, and merge_parts.sh are appropriate commands user should provide for his/her application. They are not restricted to shell codes, but can be python scripts or other executables. More specific simple examples are presented in demos/d5 and demos/d8.
The important element in this implementation is function :py:func:`.TE`.

TODO

Using external tools
====================

TODO

Working wiht GATK
=================

TODO
