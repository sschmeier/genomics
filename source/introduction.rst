Introduction
============

This is an introductory tutorial for learning genomics mostly on the Linux command-line.
Should you need to refresh your knowledge about either Linux or the command-line, have a look `here <http://linux.sschmeier.com/>`__.

In this tutorial you will learn how to analyse next-generation sequencing (NGS) data.
The data you will be using is actual research data.
The experiment follows a similar strategy as in what is called an "experimental evolution" experiment [KAWECKI2012]_, [ZEYL2006]_.
The final aim is to identify the genome variations in evolved lines of *E. coli* that can explain the observed biological phenotype(s).

The workflow
------------

The tutorial workflow is summarised in :numref:`fig-workflow`.

.. _fig-workflow:
.. figure:: /images/workflow.png

   The tutorial will follow this workflow.


Learning outcomes
-----------------

During this tutorial you will learn to:

- Check the data quality of an NGS experiment
- Create a genome assembly of the ancestor based on NGS data
- Map NGS reads of evolved lines to the created ancestral reference genome
- Call genome variations/mutations in the evolved lines
- Annotate a newly derived reference genome
- Find variants of interest that may be responsible for the observed evolved phenotype(s)

  
.. only:: html

   .. rubric:: References

               
.. [KAWECKI2012] Kawecki TJ et al. Experimental evolution. `Trends in Ecology and Evolution (2012) 27:10 <http://dx.doi.org/10.1016/j.tree.2012.06.001>`__
               
.. [ZEYL2006] Zeyl C. Experimental evolution with yeast. `FEMS Yeast Res, 2006, 685–691 <http://doi.org/10.1111/j.1567-1364.2006.00061.x>`__

