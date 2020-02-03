.. _ngs-orthology:

Orthology and Phylogeny
=======================

Preface
-------

In this section you will use some software to find orthologue genes and do phylogenetic reconstructions.


Learning outcomes
-----------------

After studying this tutorial you should be able to:

#. Use bioinformatics software to find orthologues in the NCBI database.
#. Use bioinformatics software to perform sequence alignment.
#. Use bioinformatics software to perform phylogenetic reconstructions.


Before we start
---------------

Lets see how our directory structure looks so far:

.. code:: sh

    $ cd ~/analysis
    $ ls -1F

.. code:: sh

    annotation/
    assembly/
    data/
    kraken/
    mappings/
    multiqc_data/
    trimmed/
    trimmed-fastqc/
    variants/


Make a directory for the phylogeny results (in your analysis directory):

.. code:: sh

    $ mkdir phylogeny

.. todo:: 

    Change this to gene of E.Coli

Download the fasta file of the *S. cerevisiase* TEF2 gene to the phylogeny folder:


.. code:: sh

    $ cd phylogeny
    $ curl -O http://compbio.massey.ac.nz/data/203341/s_cerev_tef2.fas


.. note::

   Should the download fail, download manually from :ref:`downloads`.



Installing the software
-----------------------


.. code:: sh

    $ conda create -n phylo blast muscle raxml-ng

This will install a |blast| executable that you can use to remotely query the NCBI database and the |muscle| alignment program that you can use to align nucleotide or protein sequences. 
We also install |raxml|, a phylogenetic tree inference tool, which uses maximum-likelihood (ML) optimality criterion. 

Finding orthologues using BLAST
-------------------------------

We will first make a |blast| database of our current assembly so that we can find the orthologous sequence of the *E.coli* gene.
To do this, we run the command ``makeblastdb``:


.. code:: sh

    # create blast db
    $ makeblastdb –in ../assembly/scaffolds.fasta –dbtype nucl


To run |blast|, we give it:

- ``-db``: The name of the database that we are BLASTing
- ``-query``: A fasta format input file
- A name for the output files
- Some notes about the format we want

First, we blast without any formatting:

.. todo:: 

    Change this to gene of E.Coli

.. code:: sh

    $ blastn –db ../assembly/scaffolds.fasta –query s_cerev_tef2.fas > blast.out


This should output a file with a set of |blast| hits similar to what you might see on the |blast| web site.

Read through the output (e.g. using ``nano``) to see what the results of your |blast| run was.

Next we will format the output a little so that it is easier to deal with.

.. todo:: 

    Change this to gene of E.Coli

.. code:: sh

    $ blastn –db ../assembly/scaffolds.fasta –query s_cerev_tef2.fas 
             –evalue 1e-100 
             –outfmt “6 length sseq” > blast_formatted.out


This will yield a file that has only the sequences of the subject, so that we can later add those to other fasta files.
However, the formatting is not perfect.
To adjust the format such that it is fasta format, open the file in an editor (e.g. ``nano``) and edit the first line so that it has a name for your sequence.
You should know the general format of a fasta-file (e.g. the first line start with a “>”).


.. hint::

   To edit in ``vi`` editor, you will need to press the escape key and “a” or “e”.
   To save in ``vi``, you will need to press the escape key and “w” (write).
   To quit ``vi``, you will need to press the escape key and “q” (quit).

Next, you have to replace the dashes (signifying indels in the |blast| result).
This can easily be done in ``vi``:
Press the escape key, followed by: ``:%s/\-//g``

Now we will |blast| a remote database to get a list of hits that are already in the NCBI database.


.. note::

   It turns out you may not be able to access this database from within the Linux distribution. In such a case, download the file named ``blast.fas`` and place it into your ``~/analysis/phylogeny/`` directory.


.. code:: sh

    $ curl -O http://compbio.massey.ac.nz/data/203341/blast_u.fas


Append the fasta file of your yeast sequence to this file, using whatever set of commands you wish/know.


.. note::

   Should the download fail, download manually from :ref:`downloads`.


Performing an alignment
-----------------------

We will use |muscle| to perform our alignment on all the sequences in the |blast| fasta file.
This syntax is very simple (change the filenames accordingly):


.. code:: sh

    $ muscle –in infile.fas –out your_alignment.aln


Building a phylogeny
--------------------

We will use |raxml| to build our phylogeny.
This uses a maximum likelihood method to infer parameters of evolution and the topology of the tree.
Again, the syntx of the command is fairly simple, except you must make sure that you are using the directory in which |raxml| sits.


The arguments are:

- ``-s``: an alignment file
- ``-m``: a model of evolution. In this case we will use a general time reversible model with gamma distributed rates (GTR+GAMMA)
- ``-n``: outfile-name
- ``-p``: specify a random number seed for the parsimony inferences

  
.. code:: sh

    $ raxmlHPC -s your_alignment.aln -m GTRGAMMA –n ecoli_tree –p 12345


Visualizing the phylogeny
-------------------------

We will use the online software `Interactive Tree of Life (iTOL) <http://itol.embl.de/upload.cgi>`__ to visualize the tree.
Navigate to this homepage.
Open the file containing your tree (``*bestTree.out``), copy the contents, and paste into the web page (in the Tree text box).

You should then be able to zoom in and out to see where your yeast taxa is.
To find out the closest relative, you will have to use the `NCBI taxa page <https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi>`__.


.. todo:: 

    Change this to E.Coli


.. todo::

   Are you certain that the yeast are related in the way that the phylogeny suggests? Why might the topology of this phylogeny not truly reflect the evolutionary history of these yeast species? 
