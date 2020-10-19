import logging
import os
import pysam
import tempfile

from . import config, constants, estimation, preprocessing, utils

logger = logging.getLogger(__name__)


def STAR_solo(
    fastqs,
    index_dir,
    out_dir,
    whitelist_path=None,
    n_threads=8,
    n_bins=50,
    temp_dir=None,
):
    """Align FASTQs with STARsolo.

    :param fastqs: list of path to FASTQs. Order matters -- STAR assumes the
                   UMI and barcode are in read 2
    :type fastqs: list
    :param index_dir: path to directory containing STAR index
    :type index_dir: str
    :param out_dir: path to directory to place STAR output
    :type out_dir: str
    :param whitelist_path: path to textfile containing barcode whitelist,
                           defaults to `None`
    :type whitelist_path: str, optional
    :param n_threads: number of threads to use, defaults to `8`
    :type n_threads: int, optional
    :param n_bins: number of bins to use when sorting BAM, defaults to `50`
    :type n_bins: int, optional
    :param temp_dir: STAR temporary directory, defaults to `None`, which
                     uses the system temporary directory
    :type temp_dir: str, optional

    :return: dictionary containing output files
    :rtype: dict
    """
    logger.info('Aligning the following FASTQs with STAR')
    for fastq in fastqs:
        logger.info((' ' * 8) + fastq)

    # out_dir must end with a directory separator
    if not out_dir.endswith(('/', '\\')):
        out_dir += os.path.sep

    command = [config.get_STAR_binary_path()] + config.STAR_SOLO_OPTIONS
    command += ['--genomeDir', index_dir]
    command += ['--readFilesIn'] + fastqs
    command += ['--runThreadN', n_threads]
    command += ['--outFileNamePrefix', out_dir]
    command += ['--soloCBwhitelist', whitelist_path or 'None']
    command += [
        '--outTmpDir',
        os.path.join(temp_dir, f'{tempfile.gettempprefix()}{next(tempfile._get_candidate_names())}')
    ]
    # TODO: currently uses hard-coded dropseq barcode and UMI positions
    command += [
        '--soloType',
        'CB_UMI_Simple',
        '--soloCBstart',
        1,
        '--soloCBlen',
        12,
        '--soloUMIstart',
        13,
        '--soloUMIlen',
        8,
    ]
    # Attempt to increase NOFILE limit if n_threads * n_bins is greater than
    # current limit
    requested = n_threads * n_bins
    current = utils.get_file_descriptor_limit()
    if requested > current:
        maximum = utils.get_max_file_descriptor_limit()

        logger.warning(
            f'Requested number of file descriptors ({requested}) exceeds '
            f'current maximum ({current}). Attempting to increase maximum '
            f'number of file descriptors to {maximum}.'
        )
        with utils.increase_file_descriptor_limit(maximum):
            # Modify n_bins if maximum > requested
            n_bins = min(n_bins, maximum // n_threads)
            if requested > maximum:
                logger.warning(
                    f'Requested number of file descriptors ({requested}) exceeds '
                    f'maximum possible defined by the OS ({maximum}). Reducing '
                    f'number of BAM bins from 50 to {n_bins}.'
                )
            command += ['--outBAMsortingBinsN', n_bins]
            utils.run_executable(command)
    else:
        command += ['--outBAMsortingBinsN', n_bins]
        utils.run_executable(command)

    solo_dir = os.path.join(out_dir, constants.STAR_SOLO_DIR)
    gene_dir = os.path.join(solo_dir, constants.STAR_GENE_DIR)
    raw_gene_dir = os.path.join(gene_dir, constants.STAR_RAW_DIR)
    filtered_gene_dir = os.path.join(gene_dir, constants.STAR_FILTERED_DIR)
    velocyto_dir = os.path.join(solo_dir, constants.STAR_VELOCYTO_DIR, constants.STAR_RAW_DIR)

    return {
        'bam': os.path.join(out_dir, constants.STAR_BAM_FILENAME),
        'gene': {
            'raw': {
                'barcodes': os.path.join(raw_gene_dir, constants.STAR_BARCODES_FILENAME),
                'features': os.path.join(raw_gene_dir, constants.STAR_FEATURES_FILENAME),
                'matrix': os.path.join(raw_gene_dir, constants.STAR_MATRIX_FILENAME),
            },
            'filtered': {
                'barcodes': os.path.join(filtered_gene_dir, constants.STAR_BARCODES_FILENAME),
                'features': os.path.join(filtered_gene_dir, constants.STAR_FEATURES_FILENAME),
                'matrix': os.path.join(filtered_gene_dir, constants.STAR_MATRIX_FILENAME),
            },
        },
        'velocyto': {
            'barcodes': os.path.join(velocyto_dir, constants.STAR_BARCODES_FILENAME),
            'features': os.path.join(velocyto_dir, constants.STAR_FEATURES_FILENAME),
            'matrix': os.path.join(velocyto_dir, constants.STAR_MATRIX_FILENAME),
        },
    }


def count(
    fastqs,
    index_dir,
    out_dir,
    use_corrected=False,
    quality=27,
    conversion='TC',
    p_group_by=None,
    pi_group_by=None,
    whitelist_path=None,
    n_threads=8,
    temp_dir=None,
    re=None,
):
    """
    """
    # Check memory.
    available_memory = utils.get_available_memory()
    if available_memory < config.RECOMMENDED_MEMORY:
        logger.warning(
            f'There is only {available_memory / (1024 ** 3):.2f} GB of free memory on the machine. '
            f'It is highly recommended to have at least {config.RECOMMENDED_MEMORY // (1024 ** 3)} GB '
            'free when running dynast. Continuing may cause dynast to crash with an out-of-memory error.'
        )

    STAR_out_dir = os.path.join(out_dir, constants.STAR_OUTPUT_DIR)
    STAR_solo_dir = os.path.join(STAR_out_dir, constants.STAR_SOLO_DIR)
    STAR_gene_dir = os.path.join(STAR_solo_dir, constants.STAR_GENE_DIR)
    STAR_raw_gene_dir = os.path.join(STAR_gene_dir, constants.STAR_RAW_DIR)
    STAR_filtered_gene_dir = os.path.join(STAR_gene_dir, constants.STAR_FILTERED_DIR)
    STAR_velocyto_dir = os.path.join(STAR_solo_dir, constants.STAR_VELOCYTO_DIR, constants.STAR_RAW_DIR)

    # Check if these files exist. If they do, we can skip alignment.
    STAR_result = {
        'bam': os.path.join(STAR_out_dir, constants.STAR_BAM_FILENAME),
        'gene': {
            'raw': {
                'barcodes': os.path.join(STAR_raw_gene_dir, constants.STAR_BARCODES_FILENAME),
                'features': os.path.join(STAR_raw_gene_dir, constants.STAR_FEATURES_FILENAME),
                'matrix': os.path.join(STAR_raw_gene_dir, constants.STAR_MATRIX_FILENAME),
            },
            'filtered': {
                'barcodes': os.path.join(STAR_filtered_gene_dir, constants.STAR_BARCODES_FILENAME),
                'features': os.path.join(STAR_filtered_gene_dir, constants.STAR_FEATURES_FILENAME),
                'matrix': os.path.join(STAR_filtered_gene_dir, constants.STAR_MATRIX_FILENAME),
            },
        },
        'velocyto': {
            'barcodes': os.path.join(STAR_velocyto_dir, constants.STAR_BARCODES_FILENAME),
            'features': os.path.join(STAR_velocyto_dir, constants.STAR_FEATURES_FILENAME),
            'matrix': os.path.join(STAR_velocyto_dir, constants.STAR_MATRIX_FILENAME),
        }
    }
    STAR_required = utils.flatten_dict_values(STAR_result)
    if not utils.all_exists(STAR_required) or re in config.RE_CHOICES[:1]:
        STAR_result = STAR_solo(
            fastqs,
            index_dir,
            STAR_out_dir,
            whitelist_path=whitelist_path,
            n_threads=n_threads,
            temp_dir=temp_dir,
        )
    else:
        logger.info(
            'Skipped STAR because alignment files already exist. '
            'Use the `--re` argument to run alignment again.'
        )

    # Check if BAM index exists and create one if it doesn't.
    bai_path = os.path.join(STAR_out_dir, constants.STAR_BAI_FILENAME)
    if not os.path.exists(bai_path) or re in config.RE_CHOICES[:1]:
        logger.info(f'Indexing {STAR_result["bam"]} with samtools')
        pysam.index(STAR_result['bam'], bai_path, '-@', str(n_threads))

    # Parse BAM and save results
    conversions_path = os.path.join(out_dir, constants.CONVERSIONS_FILENAME)
    index_path = os.path.join(out_dir, constants.INDEX_FILENAME)
    coverage_path = os.path.join(out_dir, constants.COVERAGE_FILENAME)
    if not utils.all_exists([conversions_path, index_path]) or re in config.RE_CHOICES[:2]:
        logger.info(f'Parsing read, conversion, coverage information from BAM to {conversions_path}, {coverage_path}')
        conversions_path, index_path, coverage_path = preprocessing.parse_all_reads(
            STAR_result['bam'],
            conversions_path,
            index_path,
            coverage_path,
            n_threads=n_threads,
            temp_dir=temp_dir,
        )
    else:
        logger.info(
            'Skipped read and conversion parsing from BAM because files '
            'already exist. Use the `--re` argument to parse BAM alignments again.'
        )

    # Count conversions and calculate mutation rates
    barcodes_path = os.path.join(out_dir, constants.BARCODES_FILENAME)
    genes_path = os.path.join(out_dir, constants.GENES_FILENAME)
    counts_path = os.path.join(out_dir, constants.COUNT_FILENAME)
    rates_path = os.path.join(out_dir, constants.RATES_FILENAME)
    counts_required = [
        barcodes_path,
        genes_path,
        counts_path,
        rates_path,
    ]
    df_counts = None
    if not utils.all_exists(counts_required) or re in config.RE_CHOICES[:3]:
        logger.info('Counting conversions and calculating mutation rates')
        barcodes_path, genes_path, counts_path = preprocessing.count_conversions(
            conversions_path,
            index_path,
            barcodes_path,
            genes_path,
            counts_path,
            use_corrected=use_corrected,
            quality=quality,
            n_threads=n_threads,
            temp_dir=temp_dir
        )

        df_counts = preprocessing.read_counts(counts_path)
        rates_path = preprocessing.calculate_mutation_rates(df_counts, rates_path, group_by=p_group_by)
    else:
        logger.info(
            'Skipped conversion counting and mutation rate calculation because files '
            'already exist. Use the `--re` argument to count conversions again.'
        )

    aggregates_dir = os.path.join(out_dir, constants.AGGREGATES_DIR)
    aggregates_paths = {
        conversion: os.path.join(aggregates_dir, f'{conversion}.csv')
        for conversion in preprocessing.CONVERSION_COLUMNS
    }
    if not utils.all_exists(aggregates_paths.values()) or re in config.RE_CHOICES[:4]:
        logger.info('Aggregating counts')
        os.makedirs(aggregates_dir, exist_ok=True)
        df_counts = df_counts if df_counts is not None else preprocessing.read_counts(counts_path)
        df_genes = preprocessing.read_genes(genes_path)
        aggregates_paths = preprocessing.aggregate_counts(df_counts, df_genes, aggregates_dir)
    else:
        logger.info(
            'Skipped count aggregation because files '
            'already exist. Use the `--re` argument to count conversions again.'
        )

    estimates_dir = os.path.join(out_dir, constants.ESTIMATES_DIR)
    p_e_path = os.path.join(estimates_dir, constants.P_E_FILENAME)
    p_c_path = os.path.join(estimates_dir, constants.P_C_FILENAME)
    aggregate_path = os.path.join(estimates_dir, f'{conversion}.csv')
    estimates_paths = [p_e_path, p_c_path, aggregate_path]
    df_aggregates = None
    value_columns = [conversion, conversion[0], 'count']
    if not utils.all_exists(estimates_paths) or re in config.RE_CHOICES[:5]:
        os.makedirs(estimates_dir, exist_ok=True)

        logger.info('Estimating average mismatch rate in old RNA')
        df_rates = preprocessing.read_rates(rates_path)
        p_e, p_e_path = estimation.estimate_p_e(df_rates, p_e_path, group_by=p_group_by)

        logger.info('Estimating average mismatch rate in new RNA')
        df_aggregates = df_aggregates if df_aggregates is not None else preprocessing.read_aggregates(
            aggregates_paths[conversion]
        )
        p_c, p_c_path, aggregate_path = estimation.estimate_p_c(
            df_aggregates,
            p_e,
            p_c_path,
            aggregate_path,
            group_by=p_group_by,
            value_columns=value_columns,
            n_threads=n_threads,
        )
    else:
        logger.info(
            'Skipped rate estimation because files '
            'already exist. Use the `--re` argument to calculate estimates again.'
        )

    pi_path = os.path.join(estimates_dir, constants.PI_FILENAME)
    if not utils.all_exists([pi_path]) or re in config.RE_CHOICES[:6]:
        logger.info('Estimating fraction of newly transcribed RNA')
        with open(STAR_result['gene']['filtered']['barcodes'], 'r') as f:
            barcodes = [line.strip() for line in f.readlines()]
        df_aggregates = df_aggregates if df_aggregates is not None else preprocessing.read_aggregates(
            aggregates_paths[conversion]
        )
        p_e = estimation.read_p_e(p_e_path, group_by=p_group_by)
        p_c = estimation.read_p_c(p_c_path, group_by=p_group_by)
        pi_path = estimation.estimate_pi(
            df_aggregates,
            p_e,
            p_c,
            pi_path,
            filter_dict={'barcode': barcodes},
            p_group_by=p_group_by,
            group_by=pi_group_by,
            value_columns=value_columns,
            n_threads=n_threads,
        )
    else:
        logger.info(
            'Skipped estimation of newly transcribed RNA because files '
            'already exist. Use the `--re` argument to calculate estimate again.'
        )

    adata_path = os.path.join(out_dir, constants.ADATA_FILENAME)
    logger.info('Splitting reads')
    pis = estimation.read_pi(pi_path, group_by=pi_group_by)
    adata = utils.read_STAR_count_matrix(
        STAR_result['gene']['filtered']['barcodes'],
        STAR_result['gene']['filtered']['features'],
        STAR_result['gene']['filtered']['matrix'],
    )
    adata = estimation.split_reads(adata, pis, group_by=pi_group_by)
    adata.write_h5ad(adata_path, compression='gzip')
