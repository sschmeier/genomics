.. _ngs-voi:

Variants-of-interest
====================

Preface
-------

In this section we will use our genome annotation of our reference and our genome variants in the evolved line to find variants that are interesting in terms of the observed biology.

.. NOTE::

    You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   


Overview
--------

The part of the workflow we will work on in this section can be viewed in :numref:`fig-workflow-voi`.

.. _fig-workflow-voi:
.. figure:: images/workflow.png
    
    The part of the workflow we will work on in this section marked in red.
   
     
Learning outcomes
-----------------

After studying this section of the tutorial you should be able to:

#. Identify variants of interests.
#. Understand how the variants might affect the observed biology in the evolved line.


Before we start
---------------

Lets see how our directory structure looks so far:


.. code:: bash

    cd ~/analysis
    ls -1F


.. code:: bash

    annotation/
    assembly/
    data/
    kraken/
    mappings/
    phylogeny/
    SolexaQA/
    SolexaQA++
    trimmed/
    trimmed-fastqc/
    trimmed-solexaqa/
    variants/

  
General comments for identifying variants-of-interest
-----------------------------------------------------


Things to consider when looking for variants-of-interest:

- The quality score of the variant call.
  
  * Do we call the variant with a higher then normal score?
    
- The mapping quality score.
  
  * How confident are we that the reads were mapped at the position correctly?
    
- The location of the SNP.
  
  * SNPs in larger contigs are probably more interesting than in tiny contigs.
  * Does the SNP overlap a coding region in the genome annotation?
    
- The type of SNP.

  * substitutions vs. indels 


SnpEff
------

We will be using |snpeff| to annotate our identified variants. The tool will tell us on to which genes we should focus further analyses.


Installing software
~~~~~~~~~~~~~~~~~~~
  
Tools we are going to use in this section and how to install them if you not have done it yet.


.. code:: bash

    # activate the env
    conda activate ngs
          
    # Install these tools into the conda environment
    # if not already installed
    conda install snpeff
    conda install genometools-genometools
  

Make a directory for the results (in your analysis directory) and change into
the directory:


.. code:: bash

    mkdir voi

    # change into the directory
    cd voi

         
Prepare SnpEff database
~~~~~~~~~~~~~~~~~~~~~~~

We need to create our own config-file for |snpeff|. Where is the ``snpEff.config``:


.. code:: bash

    find ~ -name snpEff.config
    /home/manager/miniconda3/envs/ngs/share/snpeff-4.3.1m-0/snpEff.config
    

This will give you the path to the ``snpEff.config``. It might be looking a bit different then the one shown here, depending on the version of |snpeff| that is installed.

Make a local copy of the ``snpEff.config`` and then edit it with an editor of your choice:


.. code:: bash

    cp /home/manager/miniconda3/envs/ngs/share/snpeff-4.3.1m-0/snpEff.config .
    nano snpEff.config

          
Make sure the data directory path in the ``snpEff.config`` looks like this:


.. code:: bash

    data.dir = ./data/

          
There is a section with databases, which starts like this:


.. code:: bash

    #-------------------------------------------------------------------------------
    # Databases & Genomes
    #
    # One entry per genome version. 
    #
    # For genome version 'ZZZ' the entries look like
    #	ZZZ.genome              : Real name for ZZZ (e.g. 'Human')
    #	ZZZ.reference           : [Optional] Comma separated list of URL to site/s Where information for building ZZZ database was extracted.
    #	ZZZ.chrName.codonTable  : [Optional] Define codon table used for chromosome 'chrName' (Default: 'codon.Standard')
    #
    #-------------------------------------------------------------------------------


Add the following two lines in the database section underneath these header lines:


.. code:: bash

    # my yeast genome
    yeastanc.genome : WildYeastAnc

          
Now, we need to create a local data folder called ``./data/yeastanc``.


.. code:: bash

    # create folders
    mkdir -p ./data/yeastanc


Copy our genome assembly to the newly created data folder.
The name needs to be ``sequences.fa`` or ``yeastanc.fa``:


.. code:: bash
    
    cp ../assembly/spades_final/scaffolds.fasta ./data/yeastanc/sequences.fa
    gzip ./data/yeastanc/sequences.fa

    
Copy our genome annotation to the data folder.
The name needs to be ``genes.gff`` (or ``genes.gtf`` for gtf-files).


.. code:: bash

    cp ../annotation/your_new_fungus.gff ./data/yeastanc/genes.gff
    gzip ./data/yeastanc/genes.gff


Now we can build a new |snpeff| database:


.. code:: bash

    snpEff build -c snpEff.config -gff3 -v yeastanc > snpEff.stdout 2> snpEff.stderr


.. note::
   Should this fail, due to gff-format of the annotation, we can try to convert the gff to gtf:


.. code:: bash

    # using genometools
    gt gff3_to_gtf ../annotation/your_new_fungus.gff -o ./data/yeastanc/genes.gtf
    gzip ./data/yeastanc/genes.gtf


Now, we can use the gtf annotation top build the database:


.. code:: bash
          
    snpEff build -c snpEff.config -gtf22 -v yeastanc > snpEff.stdout 2> snpEff.stderr


SNP annotation
~~~~~~~~~~~~~~

Now we can use our new |snpeff| database to annotate some variants, e.g.:


.. code:: bash

    snpEff -c snpEff.config yeastanc ../variants/evolved-6.freebayes.filtered.vcf.gz > evolved-6.freebayes.filtered.anno.vcf


|snpeff| adds ``ANN`` fields to the vcf-file entries that explain the effect of the variant.


.. note::

   If you are unable to do the annotation, you can download an annotated vcf-file from :ref:`downloads`.


Example
~~~~~~~

Lets look at one entry from the original vcf-file and the annotated one.
We are only interested in the 8th column, which contains information regarding the variant.
|snpeff| will add fields here :


.. code:: bash

    # evolved-6.freebayes.filtered.vcf (the original), column 8
    AB=0.5;ABP=3.0103;AC=1;AF=0.5;AN=2;AO=56;CIGAR=1X;DP=112;DPB=112;DPRA=0;EPP=3.16541;EPPR=3.16541;GTI=0;LEN=1;MEANALT=1;MQM=42;MQMR=42;NS=1;NUMALT=1;ODDS=331.872;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=2128;QR=2154;RO=56;RPL=35;RPP=10.6105;RPPR=3.63072;RPR=21;RUN=1;SAF=30;SAP=3.63072;SAR=26;SRF=31;SRP=4.40625;SRR=25;TYPE=snp

    # evolved-6.freebayes.filtered.anno.vcf, column 8
    AB=0.5;ABP=3.0103;AC=1;AF=0.5;AN=2;AO=56;CIGAR=1X;DP=112;DPB=112;DPRA=0;EPP=3.16541;EPPR=3.16541;GTI=0;LEN=1;MEANALT=1;MQM=42;MQMR=42;NS=1;NUMALT=1;ODDS=331.872;PAIRED=1;PAIREDR=1;PAO=0;PQA=0;PQR=0;PRO=0;QA=2128;QR=2154;RO=56;RPL=35;RPP=10.6105;RPPR=3.63072;RPR=21;RUN=1;SAF=30;SAP=3.63072;SAR=26;SRF=31;SRP=4.40625;SRR=25;TYPE=snp;ANN=T|missense_variant|MODERATE|CDS_NODE_40_length_1292_cov_29.5267_1_1292|GENE_CDS_NODE_40_length_1292_cov_29.5267_1_1292|transcript|TRANSCRIPT_CDS_NODE_40_length_1292_cov_29.5267_1_1292|protein_coding|1/1|c.664T>A|p.Ser222Thr|664/1292|664/1292|222/429||WARNING_TRANSCRIPT_INCOMPLETE,T|intragenic_variant|MODIFIER|GENE_NODE_40_length_1292_cov_29.5267_1_1292|GENE_NODE_40_length_1292_cov_29.5267_1_1292|gene_variant|GENE_NODE_40_length_1292_cov_29.5267_1_1292|||n.629A>T||||||  


When expecting the second entry, we find that |snpeff| added annotation information starting with ``ANN=T|missense_variant|...``.
If we look a bit more closely we find that the variant results in a amino acid change from a threonine to a serine (``c.664T>A|p.Ser222Thr``).
The codon for serine is ``TCN`` and for threonine is ``ACN``, so the variant in the first nucleotide of the codon made the amino acid change.

A quick protein |blast| of the CDS sequence where the variant was found (extracted from the ``genes.gff.gz``) shows that the closest hit is a translation elongation factor from a species called `Candida dubliniensis <https://en.wikipedia.org/wiki/Candida_dubliniensis>`_ another fungi.


.. _fig-blast-voi:
.. figure:: images/blast.png
    
    Results of a |blast| search of the CDS.
