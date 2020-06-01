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


Installing the software
-----------------------


.. code:: sh

    $ conda create -n phylo blast mafft raxml iqtree bedtools


This will install a |blast| executable that you can use to remotely query the NCBI database and the |mafft| alignment program that you can use to align sequences. 
We also install |raxml| and |iqtree|, phylogenetic tree inference tools, which use maximum-likelihood (ML) optimality criterion. 


The gene we are using
---------------------

We are using the the gene *gnd* (gluconate-6-phosphate dehydrogenase) as an example.
*gnd* is a highly polymorphic gene within E. coli populations, likely due to interstrain transfer and recombination. This may be a result of its proximity to the *rfb* region, which determines O antigen structure.

First, we are going to make a bed-file to get coordinates from file:


.. code:: sh

    # Figure out the chromosome and location of the `gnd` gene
    # there *must* be a better way to do this
    grep -B 1 ‘NODE_\|gnd’ ../annotation/*gbk

    # open a file to put the coordinates into using vi or nano
    vi gnd.bed

    # these are the coordinates from the contigs-file (Yours might be different), we copy into the vi/nano buffer
    NODE_42_length_35862_cov_7.082632	625	2031
    
    # safe the file and exit vi/nano


.. hint::

   To edit in ``vi`` editor, you will need to press the escape key and “a” or “e”.
   To save in ``vi``, you will need to press the escape key and “w” (write).
   To quit ``vi``, you will need to press the escape key and “q” (quit).


Use bedtools to extract the nucleotide sequence for the region:

.. code:: sh

    bedtools getfasta -fi ../assembly/contigs.fasta -bed gnd.bed > gnd.fasta 

Now, we have a fasta-file with exactly on genic region, the one from the *gnd* gene.


Finding orthologues using BLAST
-------------------------------

We will use a remote |blast| of ou gene against the *nr*-database:


.. code:: sh

    blastn -db nt -query gnd.fasta -remote -evalue 1e-100 -outfmt "6 qseqid sseqid sseq" > gnd_blast_hits.out

Some of the arguments explained:

- ``-db``: The name of the database that we are BLASTing against
- ``-query``: A fasta format input file
- ``-outfmt "6 qseqid sseqid sseq"``: Some notes about the format we want
- ``-evalue 1e-100``: An evalue cutoff for inclusion of results


Next, we are formating the result into fasta-format using the program ``awk``:

.. code:: sh

    awk 'BEGIN { OFS = "\n" } { print ">"$2, $3 }' gnd_blast_hits.out > gnd_blast_hits.fasta




Append the fasta file of your E. coli gene to this file, using whatever set of commands you wish/know. 
For example:

.. code:: sh

    # append the gnd gene from our E. coli!
    cat gnd.fasta >> gnd_blast_hits.fasta


Performing an alignment
-----------------------

We will use |mafft| to perform our alignment on all the sequences in the |blast| fasta file.
This syntax is very simple (change the filenames accordingly):

.. code:: sh

    $ mafft gnd_blast_hits.fasta > gnd_blast_hits.aln 


Building a phylogeny
--------------------

We will use |raxml| to build our phylogeny.
This uses a maximum likelihood method to infer parameters of evolution and the topology of the tree.
Again, the syntax of the command is fairly simple, except you must make sure that you are using the directory in which |raxml| sits.


The arguments are:

- ``-s``: an alignment file
- ``-m``: a model of evolution. In this case we will use a general time reversible model with gamma distributed rates (GTR+GAMMA)
- ``-n``: outfile-name
- ``-p``: specify a random number seed for the parsimony inferences

  
.. code:: sh

    $ raxmlHPC -s gnd_blast_hits.aln -m GTRGAMMA -n ecoli_tree -p 12345


We can also use |iqtree|, which provides more information than |raxml|.


.. code:: sh
    
    iqtree -s gnd_blast_hits.aln



Visualizing the phylogeny
-------------------------

We will use the online software `Interactive Tree of Life (iTOL) <http://itol.embl.de/upload.cgi>`__ to visualize the tree.
Navigate to this homepage.
Open the file containing your tree (``*bestTree.out``), copy the contents, and paste into the web page (in the Tree text box).

You should then be able to zoom in and out to see where your taxa is.
To find out the closest relative, you will have to use the `NCBI taxa page <https://www.ncbi.nlm.nih.gov/Taxonomy/TaxIdentifier/tax_identifier.cgi>`__.

