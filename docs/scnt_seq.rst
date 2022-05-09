scNT-seq
========
The single-cell metabolically labeled new RNA tagging sequencing (scNT-seq) was developed by [Qiu2020]_. It uses Drop-seq, which is a droplet-based scRNA-seq method [Macosko2015]_.

* Sequencing technology: Drop-seq
* Induced conversion: T>C

Alignment
^^^^^^^^^
Here, we assume the appropriate STAR index has already been built (see :ref:`ref`). A single sample will consist of a pair of FASTQs, one containing the cell barcode and UMI sequences and the other containing the biological cDNA sequences. Let's say these two FASTQs are :code:`barcode_umi.fastq.gz` and :code:`cdna.fastq.gz`.

.. code:: text

	dynast align -i path/to/STAR/index -o path/to/align/output -x dropseq cdna.fastq.gz barcode_umi.fastq.gz

This will run STAR alignment and output files to :code:`path/to/align/output`.

Consensus
^^^^^^^^^
Optionally, we can call consensus sequences for each UMI using :code:`dynast consensus`. This command requires the alignment BAM and the gene annotation GTF that was used to generate the STAR index.

.. code:: text

    dynast consensus -g path/to/GTF.gtf --barcode-tag CB --umi-tag UB path/to/align/output/Aligned.sortedByCoord.out.bam -o path/to/consensus/output

This will create a new BAM file named :code:`path/to/consensus/output/consensus.bam`, which you can then use in the next step in place of the original alignment BAM.

Quantification
^^^^^^^^^^^^^^
Finally, to quantify the number of labeled/unlabeled RNA, we run :code:`dynast count` with the appropriate alignment BAM and the gene annotation GTF that was used to generate the STAR index to :code:`-g`.

.. code:: text

	dynast count -g path/to/GTF.gtf --barcode-tag CB --umi-tag UB path/to/alignment.bam -o path/to/count/output --conversion TC

where :code:`path/to/alignment.bam` should be :code:`path/to/align/output/Aligned.sortedByCoord.out.bam` if you did not run :code:`dynast consensus`, or :code:`path/to/consensus/output/consensus.bam` if you did.

This will quantify all RNA species and write the count matrices to :code:`path/to/count/output/adata.h5ad`.
