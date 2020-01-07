.. _ngs-variantcalling:

Variant calling
===============

Preface
-------

In this section we will use our genome assembly based on the ancestor and call genetic variants in the evolved line [NIELSEN2011]_.

.. There is an accompanying lecture for this tutorial (`SNPs - GWAS - eQTLs introduction <http://dx.doi.org/10.6084/m9.figshare.1515026>`__).

.. NOTE`::

   You will encounter some **To-do** sections at times. Write the solutions and answers into a text-file.   


Overview
--------

The part of the workflow we will work on in this section can be viewed in :numref:`fig-workflow-var`.

.. _fig-workflow-var:
.. figure:: images/workflow.png

   The part of the workflow we will work on in this section marked in red.
   
     
Learning outcomes
-----------------

After studying this tutorial section you should be able to:

#. Use tools to call variants based on a reference genome.
#, Be able to describe what influences the calling of variants.


Before we start
---------------

Lets see how our directory structure looks so far:

.. code:: bash

          cd ~/analysis
          ls -1F

.. code:: bash

          assembly/
          data/
          kraken/
          mappings/
          trimmed/
          trimmed-fastqc/

   
Installing necessary software
-----------------------------
  
Tools we are going to use in this section and how to intall them if you not have done it yet.

.. code:: bash

          # activate the env
          conda activate ngs
          
          # Install these tools into the conda environment
          # if not already installed
          conda install samtools
          conda install bamtools
          conda install freebayes
          conda install bedtools
          conda install vcflib
          conda install rtg-tools
          conda install bcftools

          
Preprocessing
-------------

We first need to make an index of our reference genome as this is required by the SNP caller.
Given a scaffold/contig file in fasta-format, e.g. ``scaffolds.fasta`` which is located in the directory ``assembly/spades_final``, use |samtools| to do this:


.. code:: bash
          
          samtools faidx assembly/spades_final/scaffolds.fasta
   

Furthermore we need to pre-process our mapping files a bit further and create a bam-index file (``.bai``) for the bam-file we want to work with:


.. code:: bash
               
          bamtools index -in mappings/evolved-6.sorted.dedup.q20.bam


Lets also create a new directory for the variants:


.. code:: bash

          mkdir variants

          
Calling variants
----------------

SAMtools mpileup
~~~~~~~~~~~~~~~~

We use the sorted filtered bam-file that we produced in the mapping step before.

.. code:: bash

   # We first pile up all the reads and then call variants
   samtools mpileup -u -g -f assembly/spades_final/scaffolds.fasta mappings/evolved-6.sorted.dedup.q20.bam | bcftools call -v -m -O z -o variants/evolved-6.mpileup.vcf.gz
   
|samtools| mpileup parameter:

- ``-u``: uncompressed output
- ``-g``: generate genotype likelihoods in BCF format
- ``-f FILE``: faidx indexed reference sequence file
  
|bcftools| view parameter:

- ``-v``: output variant sites only
- ``-m``: alternative model for multiallelic and rare-variant calling
- ``-o``: output file-name
- ``-O z``: output type: 'z' compressed VCF

  
Freebayes
~~~~~~~~~

As an alternative we can do some variant calling with another tool called |freebayes|.
Given a reference genome scaffold file in fasta-format, e.g. ``scaffolds.fasta`` and the index in ``.fai`` format and a mapping file (.bam file) and a mapping index (.bai file), we can call variants with |freebayes| like so:

.. code:: bash

   # Now we call variants and pipe the results into a new file
   freebayes -f assembly/spades_final/scaffolds.fasta mappings/evolved-6.sorted.dedup.q20.bam | gzip > variants/evolved-6.freebayes.vcf.gz

         
Post-processing
---------------

Understanding the output files (.vcf)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Lets look at a vcf-file:

.. code:: bash

   # first 10 lines, which are part of the header
   zcat variants/evolved-6.mpileup.vcf.gz | head

          
.. code:: bash
   
   ##fileformat=VCFv4.2
   ##FILTER=<ID=PASS,Description="All filters passed">
   ##samtoolsVersion=1.3.1+htslib-1.3.1
   ##samtoolsCommand=samtools mpileup -g -f assembly/spades_final/scaffolds.fasta -o variants/evolved-6.mpileup.bcf mappings/evolved-6.sorted.q20.bam
   ##reference=file://assembly/spades_final/scaffolds.fasta
   ##contig=<ID=NODE_1_length_1419525_cov_15.3898,length=1419525>
   ##contig=<ID=NODE_2_length_1254443_cov_15.4779,length=1254443>
   ##contig=<ID=NODE_3_length_972329_cov_15.3966,length=972329>
   ##contig=<ID=NODE_4_length_951685_cov_15.4231,length=951685>
   ##contig=<ID=NODE_5_length_925222_cov_15.39,length=925222>
   ##contig=<ID=NODE_6_length_916533_cov_15.4426,length=916533>

Lets look at the variants:

.. code:: bash
               
   # remove header lines and look at top 4 entires
   zcat variants/evolved-6.mpileup.vcf.gz | egrep -v '##' | head -4

          
.. code:: bash
          
   #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  mappings/evolved-6.sorted.q20.bam
   NODE_1_length_1419525_cov_15.3898       24721   .       T       C       164     .       DP=12;VDB=0.205941;SGB=-0.680642;MQ0F=0;AC=2;AN=2;DP4=0,0,12,0;MQ=40     GT:PL   1/1:191,36,0
   NODE_1_length_1419525_cov_15.3898       157033  .       AAGAGAGAGAGAGAGAGAGAGAGA        AAGAGAGAGAGAGAGAGAGAGA  39.3328  .       INDEL;IDV=6;IMF=0.146341;DP=41;VDB=0.0813946;SGB=-0.616816;MQSB=1;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=13,17,3,3;MQ=42     GT:PL   0/1:75,0,255
   NODE_1_length_1419525_cov_15.3898       162469  .       T       C       19.609  .       DP=16;VDB=0.045681;SGB=-0.511536;RPB=0.032027;MQB=0.832553;BQB=0.130524;MQ0F=0;ICB=1;HOB=0.5;AC=1;AN=2;DP4=13,0,3,0;MQ=39        GT:PL   0/1:54,0,155


The fields in a vcf-file are described in he table (:numref:`table-vcf`) below:

.. _table-vcf:
.. table:: The vcf-file format fields.

   +-----+-----------+--------------------------------------------------------------------------------------+
   | Col | Field     | Description                                                                          |
   +=====+===========+======================================================================================+
   | 1   | CHROM     | Chromosome name                                                                      |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 2   | POS       | 1-based position. For an indel, this is the position preceding the indel.            |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 3   | ID        | Variant identifier. Usually the dbSNP rsID.                                          |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 4   | REF       | Reference sequence at POS involved in the variant. For a SNP, it is a single base.   |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 5   | ALT       | Comma delimited list of alternative seuqence(s).                                     |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 6   | QUAL      | Phred-scaled probability of all samples being homozygous reference.                  |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 7   | FILTER    | Semicolon delimited list of filters that the variant fails to pass.                  |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 8   | INFO      | Semicolon delimited list of variant information.                                     |
   +-----+-----------+--------------------------------------------------------------------------------------+
   | 9   | FORMAT    | Colon delimited list of the format of individual genotypes in the following fields.  |
   +-----+-----------+--------------------------------------------------------------------------------------+ 
   | 10+ | Sample(s) | Individual genotype information defined by FORMAT.                                   |
   +-----+-----------+--------------------------------------------------------------------------------------+


          
Statistics
~~~~~~~~~~

Now we can use it to do some statistics and filter our variant calls.

First, to prepare out vcf-file for querying we need to index it with ``tabix``:

.. code:: bash

   tabix -p vcf variants/evolved-6.mpileup.vcf.gz


- ``-p vcf``: input format 


We can get some quick stats with ``rtg vcfstats``:


.. code:: bash
               
   rtg vcfstats variants/evolved-6.mpileup.vcf.gz

   
Example output from ``rtg vcfstats``:


.. code::

   Location                     : variants/evolved-6.mpileup.vcf.gz
   Failed Filters               : 0
   Passed Filters               : 516
   SNPs                         : 399
   MNPs                         : 0
   Insertions                   : 104
   Deletions                    : 13
   Indels                       : 0
   Same as reference            : 0
   SNP Transitions/Transversions: 1.87 (286/153)
   Total Het/Hom ratio          : 3.20 (393/123)
   SNP Het/Hom ratio            : 8.98 (359/40)
   MNP Het/Hom ratio            : - (0/0)
   Insertion Het/Hom ratio      : 0.30 (24/80)
   Deletion Het/Hom ratio       : 3.33 (10/3)
   Indel Het/Hom ratio          : - (0/0)
   Insertion/Deletion ratio     : 8.00 (104/13)
   Indel/SNP+MNP ratio          : 0.29 (117/399)
   

   
However, we can also run |bcftools| to extract more detailed statistics about our variant calls:
   

.. code:: bash
               
   bcftools stats -F assembly/spades_final/scaffolds.fasta -s - variants/evolved-6.mpileup.vcf.gz > variants/evolved-6.mpileup.vcf.gz.stats


- ``-s -``: list of samples for sample stats, "-" to include all samples
- ``-F FILE``: faidx indexed reference sequence file to determine INDEL context

  
Now we take the stats and make some plots (e.g. :numref:`fig-vcfstats`) which are particular of interest if having multiple samples, as one can easily compare them. However, we are only working with one here:


.. code:: bash
   
   mkdir variants/plots
   plot-vcfstats -p variants/plots/ variants/evolved-6.mpileup.vcf.gz.stats

   
- ``-p``: The output files prefix, add a slash at the end to create a new directory.
   

.. _fig-vcfstats:
.. figure:: images/vcfstats.png
            
    Example of ``plot-vcfstats`` output.


Variant filtration
~~~~~~~~~~~~~~~~~~


Variant filtration is a big topic in itself [OLSEN2015]_.
There is no consens yet and research on how to best filter variants is ongoing.

We will do some simple filtration procedures here.
For one, we can filter out low quality reads.

Here, we only include variants that have quality > 30.


.. code:: bash

   # use rtg vcfffilter
   rtg vcffilter -q 30 -i variants/evolved-6.mpileup.vcf.gz -o variants/evolved-6.mpileup.q30.vcf.gz


- ``-i FILE``: input file
- ``-o FILE``: output file
- ``-q FLOAT``: minimal allowed quality in output.
  
   
or use |vcflib|:


.. code:: bash

   # or use vcflib
   zcat variants/evolved-6.mpileup.vcf.gz  | vcffilter -f "QUAL >= 30" | gzip > variants/evolved-6.mpileup.q30.vcf.gz z
      
- ``-f "QUAL >= 30"``: we only include variants that have been called with quality >= 30.


Quick stats for the filtered variants:
  
.. code:: bash 
          
   # look at stats for filtered 
   rtg vcfstats variants/evolved-6.mpileup.q30.vcf.gz


|freebayes| adds some extra information to the vcf-files it creates.
This allows for some more detailed filtering.
This strategy will NOT work on the |samtools| mpileup called variants
Here we filter, based on some recommendation form the developer of |freebayes|:


.. code:: bash

   zcat variants/evolved-6.freebayes.vcf.gz  | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" | gzip > variants/evolved-6.freebayes.filtered.vcf.gz


- ``QUAL > 1``: removes really bad sites
- ``QUAL / AO > 10``: additional contribution of each obs should be 10 log units (~ Q10 per read)
- ``SAF > 0 & SAR > 0``: reads on both strands
- ``RPR > 1 & RPL > 1``: at least two reads “balanced” to each side of the site

   
  
.. todo::
    
   Look at the statistics. One ratio that is mentioned in the statistics is transition transversion ratio (*ts/tv*).
   Explain what this ratio is and why the observed ratio makes sense. 


This strategy used here will do for our purposes.
However, several more elaborate filtering strategies have been explored, e.g. `here <https://github.com/ekg/freebayes#observation-filters-and-qualities>`__.



.. only:: html

   .. rubric:: References

.. [NIELSEN2011] Nielsen R, Paul JS, Albrechtsen A, Song YS. Genotype and SNP calling from next-generation sequencing data. `Nat Rev Genetics, 2011, 12:433-451 <http://doi.org/10.1038/nrg2986>`__

.. [OLSEN2015] Olsen ND et al. Best practices for evaluating single nucleotide variant calling methods for microbial genomics. `Front. Genet., 2015, 6:235. <https://doi.org/10.3389/fgene.2015.00235>`__
