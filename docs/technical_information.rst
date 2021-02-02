.. _technical_information:

Technical Information
=====================
This section details technical information of the quantification and statistical correction procedures of the dynast :code:`count` command. Descriptions of the :code:`ref` and :code:`align` commands are in :ref:`pipeline_usage`.

Quantification procedure
^^^^^^^^^^^^^^^^^^^^^^^^
The dynast quantification procedure consists of six steps:

1. :ref:`index`
2. :ref:`parse`
3. :ref:`snp`
4. :ref:`count`
5. :ref:`aggregate`
6. :ref:`estimate`
7. :ref:`split`

.. _index:

:code:`index`
'''''''''''''
The input BAM is indexed, if it does not exist. The BAM index will be generated at the same location and name as the input BAM, suffixed with :code:`.bai`.

.. _parse:

:code:`parse`
'''''''''''''
All files generated during this step is output to the :code:`0_parse` directory in the output directory (:code:`-o`).

1. All gene and transcript information are parsed from the gene annotation GTF (:code:`-g`) and saved as Python pickles :code:`genes.pkl.gz` and :code:`transcripts.pkl.gz`, respectively.
2. All aligned reads are parsed from the input BAM and output to :code:`conversions.csv` and :code:`no_conversions.csv`. The former contains a line for every conversion, and the latter contains a line for every read that does not have any conversions. Note that no conversion filtering (:code:`--quality`) is performed in this step. Two :code:`.idx` files are also output, corresponding to each of these CSVs, which are used downstream for fast parsing. RNA velocity types are also assigned in this step if :code:`--no-velocity` was not provided.

.. _snp:

:code:`snp`
'''''''''''
All files generated during this step is output to the :code:`0_snp` directory in the output directory (:code:`-o`). This step is skipped if :code:`--snp-threshold` is not specified.

1. Read coverage of the genome is computed by parsing all aligned reads from the input BAM and output to :code:`coverage.csv`.
2. SNPs are detected by calculating, for every genomic position, the fraction of reads with a conversion at that position over its coverage. If this fraction is greater than :code:`--snp-threshold`, then the genomic position is written to the output file :code:`snps.csv`. Any conversion with PHRED quality less than or equal to :code:`--quality` is not counted as a conversion.

.. _count:

:code:`count`
'''''''''''''
All files generated during this step is output to the :code:`1_count` directory in the output directory (:code:`-o`).

1. For every read, the numbers of each conversion (A>C, A>G, A>T, C>A, etc.) and nucleotide content (how many of A, C, G, T there are in the region that the read aligned to) are counted. Any SNPs provided with :code:`--snp-csv` or detected from the :ref:`snp` step are not counted. If both are present, the union is used. Additionally, Any conversion with PHRED quality less than or equal to :code:`--quality` is not counted as a conversion.
2. For UMI-based technologies, reads are deduplicated by the following order of priority: 1) read that aligns to the transcriptome, 2) read with the longest alignment, 3) read with least number of conversions. Reads are considered duplicates if they share the same barcode, UMI, and gene assignment. For plate-based technologies, read deduplication should have been performed in the alignment step (in the case of STAR, with the :code:`--soloUMIdedup Exact`), but in the case of multimapping reads, it becomes a bit more tricky. If a read is multimapping such that some alignments map to the transcriptome while some do not, the transcriptome alignment is taken (there can not be multiple transcriptome alignments, as this is a constraint within STAR). If none align to the transcriptome and the alignments are assigned to multiple genes, the read is dropped, as it is impossible to assign the read with confidence. If none align to the transcriptome and the alignments are assigned multiple velocity types, the velocity type is manually set to :code:`ambiguous` and the first alignment is kept. If none of these cases are true, the first alignment is kept. The final deduplicated/de-multimapped counts are output to :code:`counts.csv`.

.. _aggregate:

:code:`aggregate`
'''''''''''''''''
All files generated during this step is output to the :code:`2_aggregate` directory in the output directory (:code:`-o`). This step is skipped if :code:`--correction` is not specified.

1. Mutation rates for each base is calculated and output to :code:`rates.csv`.
2. For each cell and gene and for each conversion provided with :code:`--conversion`, the conversion counts are aggregated into a CSV file such that each row contains the following columns: cell barcode, gene, conversion count, nucleotide content of the original base (i.e. if the conversion is T>C, this would be T), and the number of reads that have this particular barcode-gene-conversion-content combination. This procedure is done for transcriptome reads, along with all velocity types as long as :code:`--no-velocity` was not specified. The resulting tables are written to :code:`transcriptome.csv`, :code:`spliced.csv`, :code:`unspliced.csv`, :code:`ambiguous.csv`, respectively.


.. _estimate:

:code:`estimate`
''''''''''''''''
All files generated during this step is output to the :code:`3_estimate` directory in the output directory (:code:`-o`). This step is skipped if :code:`--correction` is not specified.

1. The background conversion rate :math:`p_e` is estimated, if :code:`--p-e` was not provided (see :ref:`_background_estimation`). If :code:`--p-e` was provided, this value is used and estimation is skipped. :math:`p_e`s are written to :code:`p_e.csv`.
2. The induced conversion rate :math:`p_c` is estimated using an expectation maximization (EM) approach, for each conversion provided with :code:`--conversion` (see :ref:`_induced_rate_estimation`). :math:`p_c`s are written to :code:`p_c_{conversion}.csv` where :code:`{conversion}` is an underscore-delimited list of each conversion (because multiple conversions can be introduced in a single timepoint). This step is skipped for control samples with :code:`--control`.

.. _split:

:code:`split`
'''''''''''''
All files generated during this step is output to the output directory (:code:`-o`). This step is skipped if :code:`--control` is specified. All results are compiled into a single AnnData :code:`H5AD` file. The AnnData object contains the following:

* The transcriptome read counts in :code:`.X`.
* Unlabeled and labeled transcriptome read counts in :code:`.layers['X_unlabeled']` and :code:`.layers['X_labeled']`. If :code:`--correction` was specified, the corrected counts are in :code:`.layers['X_unlabeled_{conversion}_corrected']` and :code:`.layers['X_labeled_{conversion}_corrected']` where :code:`{conversion}` is an underscore-delimited list of each conversion provided with :code:`--conversion`. In addition, the actual estimated fractions of labeled RNA :math:`\pi` are in :code:`.layers['X_pi_{conversion}']`.
* [Only if :code:`--no-velocity` was not specified] Spliced, unspliced and ambiguous read counts in :code:`.layers['spliced']`, :code:`.layers['unspliced']` and :code:`.layers['ambiguous']`. If :code:`--correction` was specified, layers analogous to transcriptome read counts are added, with the exception of ambiguous read counts (i.e. no correction is ever performed on these reads).


.. _statistical_correction:

Statistical correction
^^^^^^^^^^^^^^^^^^^^^^

Overview
''''''''

.. _background_estimation:

Background estimation
'''''''''''''''''''''

.. _induced_rate_estimation:

Induced rate estimation
'''''''''''''''''''''''

Bayesian inference
''''''''''''''''''
