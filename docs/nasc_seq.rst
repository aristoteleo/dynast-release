.. _nasc_seq:

NASC-seq
=========
The new transcriptome alkylation-dependent scRNA-seq (NASC-seq) was developed by [Hendriks2019]_. It uses Smart-seq, which is a plate-based scRNA-seq method that provides great read coverage, compared to droplet-based methods [Picelli2013]_. Smart-seq experiments generate single or pairs of FASTQs for each cell sequenced, which dynast processes simultaneously.

* Sequencing technology: Smart-Seq2
* Induced conversion: T>C

Alignment
^^^^^^^^^
Here, we assume the appropriate STAR index has already been built (see :ref:`ref`). Since we have multiple sets of FASTQs, we need to prepare a FASTQ manifest CSV, instead of providing these as an argument to :code:`dynast align`. The manifest CSV contains three columns where the first column is a unique cell name/ID, the second column is the path to the first FASTQ, and the third is the path to the second FASTQ. For single-end reads, the third column can be a single :code:`-` character. Here is an example with two cells:

.. code:: text

	cell_1,path/to/R1.fastq.gz,path/to/R2.fastq.gz
	cell_2,path/to/R1.fastq.gz,-

Then, we use this manifest as the input to :code:`dynast align`.

.. code:: text

	dynast align -i path/to/STAR/index -o path/to/align/output -x smartseq manifest.csv

This will run STAR alignment and output files to :code:`path/to/align/output`.

Quantification
^^^^^^^^^^^^^^
The alignment BAM is generated at :code:`path/to/align/output/Aligned.sortedByCoord.out.bam`, which we provde as input to :code:`dynast count`. We also need to provide the gene annotation GTF that was used to generate the STAR index to :code:`-g`.

.. code:: text

	dynast count -g path/to/GTF.gtf --barcode-tag RG path/to/align/output/Aligned.sortedByCoord.out.bam -o path/to/count/output --conversion TC

This will quantify all RNA species and write the count matrices to :code:`path/to/count/output/adata.h5ad`.
