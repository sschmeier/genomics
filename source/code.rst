.. _ngs-code:

Coding solutions
================


QC
--

.. _code-fastp:

Code: fastp
~~~~~~~~~~~

.. code:: sh

    # run fastp like this on the ancestor:
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction  --cut_right --html trimmed/anc.fastp.html --json trimmed/anc.fastp.json --thread 2 -i data/anc_R1.fastq.gz -I data/anc_R2.fastq.gz -o trimmed/anc_R1.fastq.gz -O trimmed/anc_R2.fastq.gz

    # run the evolved samples through fastp
    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html trimmed/evol1.fastp.html --json trimmed/evol1.fastp.json --thread 2 -i data/evol1_R1.fastq.gz -I data/evol1_R2.fastq.gz -o trimmed/evol1_R1.fastq.gz -O trimmed/evol1_R2.fastq.gz

    fastp --detect_adapter_for_pe --overrepresentation_analysis --correction --cut_right --html trimmed/evol2.fastp.html --json trimmed/evol2.fastp.json --thread 2 -i data/evol2_R1.fastq.gz -I data/evol2_R2.fastq.gz -o trimmed/evol2_R1.fastq.gz -O trimmed/evol2_R2.fastq.gz


.. _code-qc1:

Code: FastQC
~~~~~~~~~~~~

*Create directory:*

.. code:: sh

    mkdir trimmed-fastqc


*Run FastQC:*

.. code:: sh

    fastqc -o trimmed-fastqc trimmed/*.fastq.gz
  

*Run MultiQC*

.. code:: sh

    multiqc trimmed-fastqc trimmed


*Open |multiqc| report html webpage:*

.. code:: sh

   firefox multiqc_report.html



Assembly
--------

.. _code-assembly1:

Code: SPAdes assembly (trimmed data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: sh

    spades.py -o assembly/spades-150/ --careful -1 trimmed/anc_R1.fastq.gz -2 trimmed/anc_R2.fastq.gz 

.. _code-assembly2:

Code: SPAdes assembly (original data)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: sh
    
    spades.py -o assembly/spades-original/ --careful -1 data/anc_R1.fastq.gz -2 data/anc_R2.fastq.gz


Code: Quast
~~~~~~~~~~~

.. code:: sh

    quast -o assembl/quast assembly/spades-150/scaffolds.fasta assembly/spades-original/scaffolds.fasta


Mapping
-------

.. _code-bwa1:

Code: BWA indexing
~~~~~~~~~~~~~~~~~~~~

*Index the genome assembly:*

.. code:: sh

   bwa index assembly/scaffolds.fasta


.. _code-bwa2:

Code: BWA mapping
~~~~~~~~~~~~~~~~~~~

*Run bwa mem:*

.. code:: sh

    # trimmed data
    bwa mem assembly/scaffolds.fasta trimmed/evol1_R1.fastq.gz trimmed/evol1_R2.fastq.gz > mappings/evol1.sam
    bwa mem assembly/scaffolds.fasta trimmed/evol2_R1.fastq.gz trimmed/evol2_R2.fastq.gz > mappings/evol2.sam


.. .. _code-bowtie1:

.. Code: Bowtie2 indexing
.. ~~~~~~~~~~~~~~~~~~~~~~

.. *Build the index:*

.. .. code:: sh

..    bowtie2-build assembly/scaffolds.fasta assembly/scaffolds


.. .. _code-bowtie2:

.. Code: Bowtie2 mapping
.. ~~~~~~~~~~~~~~~~~~~~~~

.. *Map to the genome. Use a max fragement length of 1000 bp:*

.. .. code:: sh

..    bowtie2 -X 1000 -x assembly/scaffolds -1 trimmed/evol1_R1.fastq.gz -2 trimmed/evol1_R2.fastq.gz -S mappings/evol1.sam
..    bowtie2 -X 1000 -x assembly/scaffolds -1 trimmed/evol2_R1.fastq.gz -2 trimmed/evol2_R2.fastq.gz -S mappings/evol2.sam


.. _code-map:

Code: Mapping post-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code:: sh

    #
    # Evol 1
    #

    # fixmate and compress to bam
    samtools sort -n -O sam mappings/evol1.sam | samtools fixmate -m -O bam - mappings/evol1.fixmate.bam
    rm mappings/evol1.sam
    # sort 
    samtools sort -O bam -o mappings/evol1.sorted.bam mappings/evol1.fixmate.bam
    rm mappings/evol1.fixmate.bam
    # mark duplicates
    samtools markdup -r -S mappings/evol1.sorted.bam mappings/evol1.sorted.dedup.bam
    rm mappings/evol1.sorted.bam
    # extract q20 mappers
    samtools view -h -b -q 20 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.dedup.q20.bam
    # extract unmapped
    samtools view -b -f 4 mappings/evol1.sorted.dedup.bam > mappings/evol1.sorted.unmapped.bam
    rm mappings/evol1.sorted.dedup.bam 
    # covert to fastq
    samtools fastq -1 mappings/evol1.sorted.unmapped.R1.fastq.gz -2 mappings/evol1.sorted.unmapped.R2.fastq.gz mappings/evol1.sorted.unmapped.bam
    # delete not needed files
    rm mappings/evol1.sorted.unmapped.bam

    #
    # Evol 2
    #

    samtools sort -n -O sam mappings/evol2.sam | samtools fixmate -m -O bam - mappings/evol2.fixmate.bam
    rm mappings/evol2.sam
    samtools sort -O bam -o mappings/evol2.sorted.bam mappings/evol2.fixmate.bam
    rm mappings/evol2.fixmate.bam
    samtools markdup -r -S mappings/evol2.sorted.bam mappings/evol2.sorted.dedup.bam
    rm mappings/evol2.sorted.bam
    samtools view -h -b -q 20 mappings/evol2.sorted.dedup.bam > mappings/evol2.sorted.dedup.q20.bam
    rm mappings/evol2.sorted.dedup.bam

.. _code-var:

Code: Variant calling
~~~~~~~~~~~~~~~~~~~~~

.. code::

    # index genome
    samtools faidx assembly/scaffolds.fasta
    mkdir variants

    #   
    # Evol 1
    #

    # index mappings
    bamtools index -in mappings/evol1.sorted.dedup.q20.bam

    # calling variants
    freebayes -p 1 -f assembly/scaffolds.fasta mappings/evol1.sorted.dedup.q20.bam > variants/evol1.freebayes.vcf
    # compress
    bgzip variants/evol1.freebayes.vcf
    # index
    $ tabix -p vcf variants/evol1.freebayes.vcf.gz

    # filtering
    zcat variants/evol1.freebayes.vcf.gz | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" | bgzip > variants/evol1.freebayes.filtered.vcf.gz
    tabix -p vcf variants/evol1.freebayes.filtered.vcf.gz

    #
    # Evol 2
    #

    # index mappings
    bamtools index -in mappings/evol2.sorted.dedup.q20.bam

    # calling variants
    freebayes -p 1 -f assembly/scaffolds.fasta mappings/evol2.sorted.dedup.q20.bam > variants/evol2.freebayes.vcf
    # compress
    bgzip variants/evol2.freebayes.vcf
    # index
    $ tabix -p vcf variants/evol2.freebayes.vcf.gz

    # filtering
    zcat variants/evol2.freebayes.vcf.gz | vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" | bgzip > variants/evol2.freebayes.filtered.vcf.gz
    tabix -p vcf variants/evol2.freebayes.filtered.vcf.gz