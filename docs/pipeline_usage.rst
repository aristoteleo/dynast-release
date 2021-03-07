.. _pipeline_usage:

Pipeline Usage
==============
This sections covers basic usage of dynast.

.. _ref:

Building the STAR index with :code:`ref`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Internally, dynast uses the STAR RNA-seq aligner to align reads to the genome [Dobin2013]_. Therefore, we must construct a STAR index to use for alignment. The :code:`dynast ref` command is a wrapper around the STAR's :code:`--runMode genomeGenerate` command, while also providing useful default parameters to limit memory usage, similar to `Cell Ranger <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger>`_. Existing STAR indices can be used interchangeably with ones generated through dynast. A genome FASTA and gene annotation GTF are required to build the STAR index.

.. code-block:: text

	usage: dynast ref [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] -i INDEX [-m MEMORY] fasta gtf

	Build a STAR index from a reference

	positional arguments:
	  fasta       Genomic FASTA file
	  gtf         Reference GTF file

	optional arguments:
	  -h, --help  Show this help message and exit
	  --tmp TMP   Override default temporary directory
	  --keep-tmp  Do not delete the tmp directory
	  --verbose   Print debugging information
	  -t THREADS  Number of threads to use (default: 8)
	  -m MEMORY   Maximum memory used, in GB (default: 16)

	required arguments:
	  -i INDEX    Path to the directory where the STAR index will be generated


Aligning FASTQs with :code:`align`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :code:`dynast align` command is a wrapper around STARsolo [Dobin2013]_. Dynast automatically formats the arguments to STARsolo to ensure the resulting alignment BAM contains information necessary for downstream processing.

Additionally, :code:`align` sets a more lenient alignment score cutoff by setting :code:`--outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3` because the reads are expected to have experimentally-induced conversions. The STARsolo defaults for both are :code:`0.66`. The :code:`--STAR-overrides` argument can be used to pass arguments directly to STAR.

:code:`dynast align` outputs a different set of BAM tags in the alignment BAM depending on the type of sequencing technology specified. These are described in the following subsections.

.. code-block:: text

	usage: dynast align [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] -i INDEX [-o OUT] -x TECHNOLOGY
	                    [--strand {forward,reverse,unstranded}] [-w WHITELIST] [--overwrite]
	                    [--STAR-overrides ARGUMENTS]
	                    fastqs [fastqs ...]

	Align FASTQs

	positional arguments:
	  fastqs                FASTQ files. If `-x smartseq`, this is a single manifest CSV file where the first
	                        column contains cell IDs and the next two columns contain paths to FASTQs (the third
	                        column may contain a dash `-` for single-end reads).

	optional arguments:
	  -h, --help            Show this help message and exit
	  --tmp TMP             Override default temporary directory
	  --keep-tmp            Do not delete the tmp directory
	  --verbose             Print debugging information
	  -t THREADS            Number of threads to use (default: 8)
	  --strand {forward,reverse,unstranded}
	                        Read strandedness. (default: `forward`)
	  -w WHITELIST          Path to file of whitelisted barcodes to correct to. If not provided, all barcodes are
	                        used.
	  --overwrite           Overwrite existing alignment files
	  --STAR-overrides ARGUMENTS
	                        Arguments to pass directly to STAR.

	required arguments:
	  -i INDEX              Path to the directory where the STAR index is located
	  -o OUT                Path to output directory (default: current directory)
	  -x TECHNOLOGY         Single-cell technology used. `dynast --list` to view all supported technologies

.. _umi_bam_tags:

UMI-based technologies
''''''''''''''''''''''
For UMI-based technologies (such as Drop-seq, 10X Chromium, scNT-seq), the following BAM tags are written to the alignment BAM.

* :code:`MD`
* :code:`CR`, :code:`CB` for raw and corrected barcodes
* :code:`UR`, :code:`UB` for raw and corrected UMIs

.. _plate_bam_tags:

Plate-based technologies
''''''''''''''''''''''''
For plate-based technologies (such as Smart-Seq), the following BAM tags are written to the alignment BAM.

* :code:`MD` BAM tag
* :code:`HI` BAM tag for paired reads
* :code:`RG` BAM tag, indicating the sample name


Quantifying counts with :code:`count`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:code:`dynast count` parses the alignment BAM and quantifies the four RNA species (unlabeled unspliced, unlabeled spliced, labeled unspliced, labeled spliced) and outputs the results as a ready-to-use `AnnData <https://anndata.readthedocs.io/en/latest/>`_ :code:`H5AD` file. In order to properly quantify the above four species, the alignment BAM must contain specific BAM tags, depending on what sequencing technology was used. If :code:`dynast align` was used to generate the alignment BAM, dynast automatically configures the appropriate BAM tags to be written.

.. code-block:: text

	usage: dynast count [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] -g GTF --conversion CONVERSION [-o OUT]
	                    [--umi-tag TAG] [--barcode-tag TAG] [--gene-tag TAG] [--strand {forward,reverse,unstranded}]
	                    [--quality QUALITY] [--snp-threshold THRESHOLD] [--snp-csv CSV] [--barcodes BARCODES]
	                    [--no-splicing] [--control] [--overwrite]
	                    bam

	Quantify unlabeled and labeled RNA

	positional arguments:
	  bam                   Alignment BAM file that contains the appropriate UMI and barcode tags, specifiable with
	                        `--umi-tag`, and `--barcode-tag`.

	optional arguments:
	  -h, --help            Show this help message and exit
	  --tmp TMP             Override default temporary directory
	  --keep-tmp            Do not delete the tmp directory
	  --verbose             Print debugging information
	  -t THREADS            Number of threads to use (default: 8)
	  -o OUT                Path to output directory (default: current directory)
	  --umi-tag TAG         BAM tag to use as unique molecular identifiers (UMI). If not provided, all reads are assumed
	                        to be unique. (default: None)
	  --barcode-tag TAG     BAM tag to use as cell barcodes. If not provided, all reads are assumed to be from a single
	                        cell. (default: None)
	  --gene-tag TAG        BAM tag to use as gene assignments (default: GX)
	  --strand {forward,reverse,unstranded}
	                        Read strandedness. (default: `forward`)
	  --quality QUALITY     Base quality threshold. Only bases with PHRED quality greater than this value will be
	                        considered when counting conversions. (default: 27)
	  --snp-threshold THRESHOLD
	                        Conversions with (# conversions) / (# reads) greater than this threshold will be considered a
	                        SNP and ignored. (default: no SNP detection)
	  --snp-csv CSV         CSV file of two columns: contig (i.e. chromosome) and genome position of known SNPs
	  --barcodes BARCODES   Textfile containing filtered cell barcodes. Only these barcodes will be processed.
	  --no-splicing, --transcriptome-only
	                        Do not assign reads a splicing status (spliced, unspliced, ambiguous) and ignore reads that
	                        are not assigned to the transcriptome.
	  --control             Indicate this is a control sample, which is used to detect SNPs.
	  --overwrite

	required arguments:
	  -g GTF                Path to GTF file used to generate the STAR index
	  --conversion CONVERSION
	                        The type of conversion(s) introduced at a single timepoint. Multiple conversions can be
	                        specified with a comma-delimited list. For example, T>C and A>G is TC,AG. This option can be
	                        specified multiple times (i.e. dual labeling), for each labeling timepoint.

Basic arguments
'''''''''''''''
The :code:`--barcode-tag` and :code:`--umi-tag` arguments are used to specify what BAM tags should be used to differentiate cells (barcode) and RNA molecules (UMI). If the former is not specified, all BAM alignments are assumed to be from a single cell, and if the latter is not specified, all aligned reads are assumed to be unique (i.e. no read deduplication is performed). If :code:`align` was used to generate the alignment BAM, then :code:`--barcode-tag CB --umi-tag UB` is recommended for UMI-based technologies (see :ref:`umi_bam_tags`), and :code:`--barcode-tag RG` is recommended for Plate-based technologies (see :ref:`plate_bam_tags`).

The :code:`--strand` argument can be used to specify the read strand of the sequencing technology. Usually, the default (:code:`forward`) is appropriate, but this argument may be of use for other technologies.

The :code:`--conversion` argument is used to specify the type of conversion that is experimentally introduced as a two-character string. For instance, a T>C conversion is represented as :code:`TC`, which is the default. Multiple conversions can be specified as a comma-delimited list, and :code:`--conversion` may be specified multiple times to indicate multiple-indexing experiments. For example, for an experiment that introduced T>C mutations at timepoint 1 and A>G and C>G mutations at timepoint 2, the appropriate options would be :code:`--conversion TC --conversion AG,CG`.

.. _snps:

Detecting and filtering SNPs
''''''''''''''''''''''''''''
:code:`dynast count` has the ability to detect single-nucleotide polymorphisms (SNPs) by calculating the fraction of reads with a mutation at a certain genomic position. :code:`--snp-threshold` can be used to specify the proportion threshold greater than which a SNP will be called at that position. All conversions/mutations at the genomic positions with SNPs detected in this manner will be filtered out from further processing. In addition, a CSV file containing known SNP positions can be provided with the :code:`--snp-csv` argument. This argument accepts a CSV file containing two columns: contig (i.e. chromosome) and genomic position of known SNPs.

.. _estimate:

Estimating counts with :code:`estimate`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The fraction of labeled RNA is estimated with the :code:`dynast estimate` command. Whereas :code:`dynast count` produces naive UMI count matrices, :code:`dynast estimate` statistically models labeling dynamics to estimate the true fraction of labeled RNA (and then in turn uses this fraction to split the total UMI counts into unlabeled and labeled RNA). See :ref:`statistical_estimation` of a technical overview of this process. In this section, we will simply be describing the command-line usage of this command.

.. code-block:: text

	usage: dynast estimate [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS]
	                       [--reads {total,transcriptome,spliced,unspliced}] [-o OUT] [--groups CSV]
	                       [--cell-threshold COUNT] [--cell-gene-threshold COUNT] [--control] [--p-e P_E]
	                       count_dir

	Estimate fraction of labeled RNA

	positional arguments:
	  count_dir             Path to directory that contains `dynast count` output.

	optional arguments:
	  -h, --help            Show this help message and exit
	  --tmp TMP             Override default temporary directory
	  --keep-tmp            Do not delete the tmp directory
	  --verbose             Print debugging information
	  -t THREADS            Number of threads to use (default: 8)
	  --reads {total,transcriptome,spliced,unspliced}
	                        Read groups to perform estimation on. This option can be used multiple times to estimate
	                        multiple groups. By default, all possible reads groups are used.
	  -o OUT                Path to output directory (default: current directory)
	  --groups CSV          CSV containing cell (barcode) groups, where the first column is the barcode and the second is
	                        the group name the cell belongs to. Cells will be combined per group for estimation.
	  --cell-threshold COUNT
	                        A cell must have at least this many reads for estimation. (default: 1000)
	  --cell-gene-threshold COUNT
	                        A cell-gene pair must have at least this many reads for estimation. (default: 16)
	  --control             Indicate this is a control sample, only the background mutation rate will be estimated.
	  --p-e P_E             Textfile containing a single number, indicating the estimated background mutation rate

Estimation thresholds
'''''''''''''''''''''
The :code:`--cell-threshold` and :code:`--cell-gene-threshold` arguments control the minimum number of reads that a cell and cell-gene combination must have for accurate estimation. By default, these are :code:`1000` and :code:`16` respectively. Any cells with reads less than the former are excluded from estimation, and the same goes for any genes within a cell that has less reads than the latter. If :code:`--groups` is also provided, then these thresholds apply to each cell **group** instead of each cell individually. Internally, :code:`--cell-threshold` is used to filter cells before estimating the average conversion rate in labeled RNA (see :ref:`induced_rate_estimation`), and :code:`--cell-gene-threshold` is used to filter cell-gene combinations before estimating the fraction of new RNA (see :ref:`bayesian_inference`).

Estimation on a subset of RNA species
'''''''''''''''''''''''''''''''''''''
The :code:`--reads` argument controls which RNA species to run the estimation procedure on. By default, all possible RNA species, minus :code:`ambiguous` reads, are used. This argument can take on the following values: :code:`total`, :code:`transcriptome`, :code:`spliced`, :code:`unspliced` (see :ref:`read_groups`). The value of this argument specifies which group of unlabeled/labeled RNA counts will be estimated. For instance, :code:`--reads spliced` will run statistical estimation on unlabeled/labeled spliced reads. This option may be provided multiple times to run estimation on multiple groups. The procedure involves estimating the conversion rate of unlabeled and labeled RNA, and modeling the fraction of new RNA as a binomial mixture model (see :ref:`statistical_estimation`).


Control samples
^^^^^^^^^^^^^^^
Control samples may be used to find common SNPs and directly estimate the conversion rate of unlabeled RNA (see :ref:`background_estimation`). Normally, the latter is estimating using the reads directly. However, it is possible to use a control sample (prepared in absence of the experimental introduction of conversions) to calculate this value directly. In addition, SNPs can be called in the control sample, and these called SNPs can be used when running the test sample(s) (see :ref:`snps` for SNP arguments). Note that SNP calling is done with :code:`dynast count`.

A typical workflow for a control sample is the following.

.. code-block:: text

  dynast count --control --snp-threshold 0.5 -o control_count --conversion TC -g GTF.gtf CONTROL.bam
	dynast estimate --control -o control_estimate control_count

The :code:`dynast count` command detects SNPs from the control sample and outputs them to the file :code:`snps.csv` in the output directory :code:`control_count`. The :code:`dynast estimate` calculates the background conversion rate of unlabeled RNA to the file :code:`p_e.csv` in the output directory :code:`control_estimate`. These files can then be used as input when running the test sample.

.. code-block:: text

	dynast count --snp-csv control_count/snps.csv -o test_count --conversion TC -g GTF.gtf INPUT.bam
	dynast estimate --p-e control_estimate/p_e.csv -o test_estimate test_count

The above set of commands runs quantification and estimation on the test sample using the SNPs detected from the control sample (:code:`control_count/snps.csv`) and the background conversion rate estimated from the control sample (:code:`control_estimate/p_e.csv`).
