sci-fate
========
The single-cell combinatorial indexing and messenger RNA labeling (sci-fate) was developed by [Cao2020]_.

* Sequencing technology: sci-fate
* Induced conversion: T>C

Alignment
^^^^^^^^^
Here, we assume the appropriate STAR index has already been built (see :ref:`ref`). A single sample will consist of a pair of FASTQs, one containing the cell barcode and UMI sequences and the other containing the biological cDNA sequences. Let's say these two FASTQs are :code:`barcode_umi.fastq.gz` and :code:`cdna.fastq.gz`.

.. code:: text

	dynast align -i path/to/STAR/index -o path/to/align/output -x scifate cdna.fastq.gz barcode_umi.fastq.gz

This will run STAR alignment and output files to :code:`path/to/align/output`.

Quantification
^^^^^^^^^^^^^^
The alignment BAM is generated at :code:`path/to/align/output/Aligned.sortedByCoord.out.bam`, which we provde as input to :code:`dynast count`. We also need to provide the gene annotation GTF that was used to generate the STAR index to :code:`-g`.

.. code:: text

	dynast count -g path/to/GTF.gtf --barcode-tag CB --umi-tag UB path/to/align/output/Aligned.sortedByCoord.out.bam -o path/to/count/output --conversion TC

This will quantify all RNA species and write the count matrices to :code:`path/to/count/output/adata.h5ad`.

.. [Cao2020] https://doi.org/10.1038/s41587-020-0480-9
