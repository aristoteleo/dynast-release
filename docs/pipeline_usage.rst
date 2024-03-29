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
* :code:`HI`, :code:`AS` for alignment index and score
* :code:`CR`, :code:`CB` for raw and corrected barcodes
* :code:`UR`, :code:`UB` for raw and corrected UMIs

.. _plate_bam_tags:

Plate-based technologies
''''''''''''''''''''''''
For plate-based technologies (such as Smart-Seq), the following BAM tags are written to the alignment BAM. See

* :code:`MD`
* :code:`HI`, :code:`AS` for alignment index and score
* :code:`RG` indicating the sample name

Calling consensus sequences with :code:`consensus`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:code:`dynast consensus` parses the alignment BAM to generate consensus sequences for each sequenced mRNA molecule (see :ref:`consensus_procedure`).

.. code-block:: text

    usage: dynast consensus [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] -g GTF [-o OUT] [--umi-tag TAG]
                            [--barcode-tag TAG] [--gene-tag TAG] [--strand {forward,reverse,unstranded}]
                            [--quality QUALITY] [--barcodes TXT] [--add-RS-RI]
                            bam

    Generate consensus sequences

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
      --quality QUALITY     Base quality threshold. When generating a consensus nucleotide at a certain position, the base
                            with smallest error probability below this quality threshold is chosen. If no base meets this
                            criteria, the reference base is chosen. (default: 27)
      --barcodes TXT        Textfile containing filtered cell barcodes. Only these barcodes will be processed.
      --add-RS-RI           Add custom RS and RI tags to the output BAM, each of which contain a semi-colon delimited list
                            of read names (RS) and alignment indices (RI) of the reads and alignments from which the
                            consensus is derived. This option is useful for debugging.

    required arguments:
      -g GTF                Path to GTF file used to generate the STAR index

The resulting BAM will contain a collection of consensus alignments and a subset of original alignments (for those alignments for which a consensus could not be determined). For the latter, they will contain the following modified BAM tags.

* :code:`GX`, :code:`GN` each containing the assigned gene ID and name. Note that these tags are used regardless of what was provided to :code:`--gene-tag`. Since these are reads for which a consensus could not be determined, these tags will be identical to what was contained in the tag provided with :code:`--gene-tag`.

For the consensus reads, their names will be seemingly random sequences of letters and numbers (in reality, these are SHA256 checksums of the grouped read names). They will also contain the following modified BAM tags.

* :code:`AS` is now the *sum* of the alignment scores of the reads
* :code:`HI`, the alignment index, is always 1

and the follwing additional BAM tags.

* :code:`RN` indicating how many reads were used to generate the consensus
* :code:`RS`, :code:`RI` each containing a semicolon-delimited list of read names and their corresponding alignment indices (:code:`HI` tag in the original BAM) that were used to generate the consensus (only added if :code:`--add-RS-RI` is provided)
* :code:`GX`, :code:`GN` each containing the assigned gene ID and name. Note that these tags are used regardless of what was provided to :code:`--gene-tag`.

Quantifying counts with :code:`count`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:code:`dynast count` parses the alignment BAM and quantifies the four RNA species (unlabeled unspliced, unlabeled spliced, labeled unspliced, labeled spliced) and outputs the results as a ready-to-use `AnnData <https://anndata.readthedocs.io/en/latest/>`_ :code:`H5AD` file. In order to properly quantify the above four species, the alignment BAM must contain specific BAM tags, depending on what sequencing technology was used. If :code:`dynast align` was used to generate the alignment BAM, dynast automatically configures the appropriate BAM tags to be written.

.. code-block:: text

    usage: dynast count [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] -g GTF --conversion CONVERSION [-o OUT]
                    [--umi-tag TAG] [--barcode-tag TAG] [--gene-tag TAG] [--strand {forward,reverse,unstranded,auto}]
                    [--quality QUALITY] [--snp-threshold THRESHOLD] [--snp-min-coverage THRESHOLD] [--snp-csv CSV]
                    [--barcodes TXT] [--gene-names] [--no-splicing | --exon-overlap {lenient,strict}] [--control]
                    [--dedup-mode {auto,conversion,exon}] [--overwrite]
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
      --strand {forward,reverse,unstranded,auto}
                            Read strandedness. By default, this is auto-detected from the BAM.
      --quality QUALITY     Base quality threshold. Only bases with PHRED quality greater than this value will be
                            considered when counting conversions. (default: 27)
      --snp-threshold THRESHOLD
                            Conversions with (# conversions) / (# reads) greater than this threshold will be considered a
                            SNP and ignored. (default: no SNP detection)
      --snp-min-coverage THRESHOLD
                            For a conversion to be considered as a SNP, there must be at least this many reads mapping to
                            that region. (default: 1)
      --snp-csv CSV         CSV file of two columns: contig (i.e. chromosome) and genome position of known SNPs
      --barcodes TXT        Textfile containing filtered cell barcodes. Only these barcodes will be processed.
      --gene-names          Group counts by gene names instead of gene IDs when generating the h5ad file.
      --no-splicing, --transcriptome-only
                            Do not assign reads a splicing status (spliced, unspliced, ambiguous) and ignore reads that
                            are not assigned to the transcriptome.
      --exon-overlap {lenient,strict}
                            Algorithm to use to detect spliced reads (that overlap exons). May be `strict`, which assigns
                            reads as spliced if it only overlaps exons, or `lenient`, which assigns reads as spliced if it
                            does not overlap with any introns of at least one transcript. (default: strict)
      --control             Indicate this is a control sample, which is used to detect SNPs.
      --dedup-mode {auto,conversion,exon}
                            Deduplication mode for UMI-based technologies (required `--umi-tag`). Available choices are:
                            `auto`, `conversion`, `exon`. When `conversion` is used, reads that have at least one of the
                            provided conversions is prioritized. When `exon` is used, exonic reads are prioritized. By
                            default (`auto`), the BAM is inspected to select the appropriate mode.
      --overwrite           Overwrite existing files.

    required arguments:
      -g GTF                Path to GTF file used to generate the STAR index
      --conversion CONVERSION
                            The type of conversion(s) introduced at a single timepoint. Multiple conversions can be
                            specified with a comma-delimited list. For example, T>C and A>G is TC,AG. This option can be
                            specified multiple times (i.e. dual labeling), for each labeling timepoint.

.. _basic_arguments:

Basic arguments
'''''''''''''''
The :code:`--barcode-tag` and :code:`--umi-tag` arguments are used to specify what BAM tags should be used to differentiate cells (barcode) and RNA molecules (UMI). If the former is not specified, all BAM alignments are assumed to be from a single cell, and if the latter is not specified, all aligned reads are assumed to be unique (i.e. no read deduplication is performed). If :code:`align` was used to generate the alignment BAM, then :code:`--barcode-tag CB --umi-tag UB` is recommended for UMI-based technologies (see :ref:`umi_bam_tags`), and :code:`--barcode-tag RG` is recommended for Plate-based technologies (see :ref:`plate_bam_tags`).

The :code:`--strand` argument can be used to specify the read strand of the sequencing technology. Usually, the default (:code:`forward`) is appropriate, but this argument may be of use for other technologies.

The :code:`--conversion` argument is used to specify the type of conversion that is experimentally introduced as a two-character string. For instance, a T>C conversion is represented as :code:`TC`, which is the default. Multiple conversions can be specified as a comma-delimited list, and :code:`--conversion` may be specified multiple times to indicate multiple-indexing experiments. For example, for an experiment that introduced T>C mutations at timepoint 1 and A>G and C>G mutations at timepoint 2, the appropriate options would be :code:`--conversion TC --conversion AG,CG`.

The :code:`--gene-names` argument can be used to specify that the resulting AnnData should contain gene names as its columns, instead of the usual gene IDs.

.. _snps:

Detecting and filtering SNPs
''''''''''''''''''''''''''''
:code:`dynast count` has the ability to detect single-nucleotide polymorphisms (SNPs) by calculating the fraction of reads with a mutation at a certain genomic position. :code:`--snp-threshold` can be used to specify the proportion threshold greater than which a SNP will be called at that position. All conversions/mutations at the genomic positions with SNPs detected in this manner will be filtered out from further processing. In addition, a CSV file containing known SNP positions can be provided with the :code:`--snp-csv` argument. This argument accepts a CSV file containing two columns: contig (i.e. chromosome) and genomic position of known SNPs.

Read deduplication modes
''''''''''''''''''''''''
The :code:`--dedup-mode` option is used to select how duplicate reads should be deduplicated for UMI-based technologies (i.e. :code:`--umi-tag` is provided). Two different modes are supported: :code:`conversion` and :code:`exon`. The former prioritizes reads that have at least one conversions provided by :code:`--conversion`. The latter prioritizes exonic reads. See :ref:`quant` for a more technical description of how deduplication is performed. Additionally, see :ref:`consensus_procedure` to get an idea of why selecting the correct option may be important.

By default, the :code:`--dedup-mode` is set to :code:`auto`, which sets the deduplication mode to :code:`exon` if the input BAM is detected to be a consensus-called BAM (a BAM generated with :code:`dynast consensus`). Otherwise, it is set to :code:`conversion`. This option has no effect for non-UMI technologies.

.. _estimate:

Estimating counts with :code:`estimate`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The fraction of labeled RNA is estimated with the :code:`dynast estimate` command. Whereas :code:`dynast count` produces naive UMI count matrices, :code:`dynast estimate` statistically models labeling dynamics to estimate the true fraction of labeled RNA (and then in turn uses this fraction to split the total UMI counts into unlabeled and labeled RNA). See :ref:`statistical_estimation` of a technical overview of this process. In this section, we will simply be describing the command-line usage of this command.

.. code-block:: text

    usage: dynast estimate [-h] [--tmp TMP] [--keep-tmp] [--verbose] [-t THREADS] [--reads {total,transcriptome,spliced,unspliced}] [-o OUT] [--barcodes TXT]
                       [--groups CSV] [--method {pi_g,alpha}] [--ignore-groups-for-est] [--genes TXT] [--cell-threshold COUNT] [--cell-gene-threshold COUNT]
                       [--gene-names] [--downsample NUM] [--downsample-mode MODE] [--control] [--p-e P_E]
                       count_dirs [count_dirs ...]

    Estimate fraction of labeled RNA

    positional arguments:
      count_dirs            Path to directory that contains `dynast count` output. When multiple are provided, the barcodes in each of the count directories are
                            suffixed with `-i` where i is a 0-indexed integer.

    optional arguments:
      -h, --help            Show this help message and exit
      --tmp TMP             Override default temporary directory
      --keep-tmp            Do not delete the tmp directory
      --verbose             Print debugging information
      -t THREADS            Number of threads to use (default: 8)
      --reads {total,transcriptome,spliced,unspliced}
                            Read groups to perform estimation on. This option can be used multiple times to estimate multiple groups. (default: all possible reads
                            groups)
      -o OUT                Path to output directory (default: current directory)
      --barcodes TXT        Textfile containing filtered cell barcodes. Only these barcodes will be processed. This option may be used multiple times when
                            multiple input directories are provided.
      --groups CSV          CSV containing cell (barcode) groups, where the first column is the barcode and the second is the group name the cell belongs to.
                            Estimation will be performed by aggregating UMIs per group. This option may be used multiple times when multiple input directories are
                            provided.
      --method {pi_g,alpha}
                            Correction method to use. May be `pi_g` to estimate the fraction of labeled RNA for every cell-gene combination, or `alpha` to use
                            alpha correction as used in the scNT-seq paper. `alpha` is recommended for UMI-based assays. This option has no effect when used with
                            `--control`. (default: alpha)
      --ignore-groups-for-est
                            Ignore cell groupings when calculating final estimations for the fraction of labeled RNA. When `--method pi_g`, groups are ignored
                            when estimating fraction of labeled RNA. When `--method alpha`, groups are ignored when estimating detection rate. This option only
                            has an effect when `--groups` is also specified.
      --genes TXT           Textfile containing list of genes to use. All other genes will be treated as if they do not exist.
      --cell-threshold COUNT
                            A cell must have at least this many reads for correction. (default: 1000)
      --cell-gene-threshold COUNT
                            A cell-gene pair must have at least this many reads for correction. Only for `--method pi_g`. (default: 16)
      --gene-names          Group counts by gene names instead of gene IDs when generating H5AD file
      --downsample NUM      Downsample the number of reads (UMIs). If a decimal between 0 and 1 is given, then the number is interpreted as the proportion of
                            remaining reads. If an integer is given, the number is interpreted as the absolute number of remaining reads.
      --downsample-mode MODE
                            Downsampling mode. Can be one of: `uniform`, `cell`, `group`. If `uniform`, all reads (UMIs) are downsampled uniformly at random. If
                            `cell`, only cells that have more reads than the argument to `--downsample` are downsampled to exactly that number. If `group`,
                            identical to `cell` but per group specified by `--groups`.
      --control             Indicate this is a control sample, only the background mutation rate will be estimated.
      --p-e P_E             Textfile containing a single number, indicating the estimated background mutation rate

Estimation methods
''''''''''''''''''
Dynast supports two different statistical correction methods. The :code:`--method pi_g` employs a Bayesian inference approach to directly estimate the fraction of labeled RNA for each cell-gene combination. While this approach performs well for plate-based assays (such as those using Smart-Seq), droplet-based assays (such as those using Drop-seq) produce very sparse counts for which this estimation procedure often fails due to low number of reads per cell-gene. Therefore, :code:`--method alpha` uses the detection rate estimation used in [Qiu2020]_, which is more suited for sparse data. See :ref:`bayesian_inference` for more information.

Estimation thresholds
'''''''''''''''''''''
The :code:`--cell-threshold` and :code:`--cell-gene-threshold` arguments control the minimum number of reads that a cell and cell-gene combination must have for accurate estimation. By default, these are :code:`1000` and :code:`16` respectively. Any cells with reads less than the former are excluded from estimation, and the same goes for any genes within a cell that has less reads than the latter. If :code:`--groups` is also provided, then these thresholds apply to each cell **group** instead of each cell individually. Internally, :code:`--cell-threshold` is used to filter cells before estimating the average conversion rate in labeled RNA (see :ref:`induced_rate_estimation`), and :code:`--cell-gene-threshold` is used to filter cell-gene combinations before estimating the fraction of new RNA and only has an effect when :code:`--method pi_g` (see :ref:`bayesian_inference`).

Estimation on a subset of RNA species
'''''''''''''''''''''''''''''''''''''
The :code:`--reads` argument controls which RNA species to run the estimation procedure on. By default, all possible RNA species, minus :code:`ambiguous` reads, are used. This argument can take on the following values: :code:`total`, :code:`transcriptome`, :code:`spliced`, :code:`unspliced` (see :ref:`read_groups`). The value of this argument specifies which group of unlabeled/labeled RNA counts will be estimated. For instance, :code:`--reads spliced` will run statistical estimation on unlabeled/labeled spliced reads. This option may be provided multiple times to run estimation on multiple groups. The procedure involves estimating the conversion rate of unlabeled and labeled RNA, and modeling the fraction of new RNA as a binomial mixture model (see :ref:`statistical_estimation`).

Grouping cells
''''''''''''''
Sometimes, grouping read counts across cells may provide better estimation results, especially in the case of droplet-based methods, which result in fewer reads per cell and gene compared to plate-based methods. The :code:`--groups` argument can be used to provide a CSV of two columns: the first containing the cell barcodes and the second containing group names that each cell belongs to. Estimation is then performed on a per-group basis by combining the read counts across all cells in each group. This strategy may be applied across different samples, simply by specifying multiple input directories. In this case, the number of group CSVs specified with :code:`--groups` must match the number of input directories. For example, when providing two input directories :code:`./input1` and :code:`./input2`, with the intention of grouping cells across these two samples, two group CSVs, :code:`groups1.csv` and :code:`groups2.csv` must be provided where the former are groups for barcodes in the first sample, and the latter are groups for barcodes in the second sample. The group names may be shared across samples. The output AnnData will still contain reads per cell.

Cell groupings provided this way may be ignored for estimation of the fraction of labeled RNA when :code:`--method pi_g` or the detection rate when :code:`--method alpha` (see :ref:`bayesian_inference`) by providing the :code:`--ignore-groups-for-est` flag. This flag may be used only in conjunction with :code:`--groups`, and when it is provided, final estimation is performed per cell, while estimation of background and induced mutation rates are still done per group.

Downsampling
''''''''''''
Downsampling UMIs uniformly, per cell, or per cell group may be useful to significantly reduce runtime while troubleshooting pipeline parameters (or just to quickly get some preliminary results). Dynast can perform downsampling when the :code:`--downsample` argument is used. The value of this argument may either be an integer indicating the number of UMIs to retain or a proportion between 0 and 1 indicating the proportion of UMIs to retain. Additionally, the downsampling mode may be specified with the :code:`--downsample-mode` argument, which takes one of the following three parameters: :code:`uniform`, :code:`cell`, :code:`group`. :code:`uniform` is the default that downsamples UMIs uniformly at random. When :code:`cell` is provided, the value of :code:`--downsample` may only be an integer specifying the threshold to downsample cells to. Only cells with UMI counts greater than this value will be downsampled to exactly this value. :code:`group` works the same way, but for cell groups and may be used only in conjunction with :code:`--groups`.

Control samples
^^^^^^^^^^^^^^^
Control samples may be used to find common SNPs and directly estimate the conversion rate of unlabeled RNA (see :ref:`background_estimation`). Normally, the latter is estimating using the reads directly. However, it is possible to use a control sample (prepared in absence of the experimental introduction of conversions) to calculate this value directly. In addition, SNPs can be called in the control sample, and these called SNPs can be used when running the test sample(s) (see :ref:`snps` for SNP arguments). Note that SNP calling is done with :code:`dynast count`.

A typical workflow for a control sample is the following.

.. code-block:: text

    dynast count --control --snp-threshold 0.5 [...] -o control_count --conversion TC -g GTF.gtf CONTROL.bam
    dynast estimate --control -o control_estimate control_count

Where :code:`[...]` indicates the usual options that would be used for :code:`dynast count` if this were not control samples. See :ref:`basic_arguments` for these options.

The :code:`dynast count` command detects SNPs from the control sample and outputs them to the file :code:`snps.csv` in the output directory :code:`control_count`. The :code:`dynast estimate` calculates the background conversion rate of unlabeled RNA to the file :code:`p_e.csv` in the output directory :code:`control_estimate`. These files can then be used as input when running the test sample.

.. code-block:: text

    dynast count --snp-csv control_count/snps.csv -o test_count [...] INPUT.bam
    dynast estimate --p-e control_estimate/p_e.csv -o test_estimate test_count

The above set of commands runs quantification and estimation on the test sample using the SNPs detected from the control sample (:code:`control_count/snps.csv`) and the background conversion rate estimated from the control sample (:code:`control_estimate/p_e.csv`).
