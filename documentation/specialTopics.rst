**************
Special topics
**************

Wokring with clusters
=====================

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

Explanation.

TODO

Using external tools
====================

TODO

Working wiht GATK
=================

TODO
