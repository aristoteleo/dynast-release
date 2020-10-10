import logging
import os
import pysam
import tempfile

from . import bam, config, constants, conversions, utils

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

    :return:
    :rtype:
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
    quality=27,
    whitelist_path=None,
    n_threads=8,
    temp_dir=None,
    realign=False,
    reparse=False,
    recount=False,
):
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
    found = [path for path in STAR_required if os.path.exists(path)]
    not_found = [path for path in STAR_required if not os.path.exists(path)]
    logger.debug(f'found: {found}')
    logger.debug(f'not found: {not_found}')
    if not_found or realign:
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
            'Skipping STAR because alignment files already exist. '
            'Use the `--realign` flag to run alignment again.'
        )

    # Check if BAM index exists and create one if it doesn't.
    bai_path = os.path.join(STAR_out_dir, constants.STAR_BAI_FILENAME)
    if not os.path.exists(bai_path) or realign:
        logger.info(f'Indexing {STAR_result["bam"]} with samtools')
        pysam.index(STAR_result['bam'], bai_path, '-@', str(n_threads))

    # Parse BAM and save results
    conversions_path = os.path.join(out_dir, constants.CONVERSIONS_FILENAME)
    index_path = os.path.join(out_dir, constants.INDEX_FILENAME)
    coverage_path = os.path.join(out_dir, constants.COVERAGE_FILENAME)
    if not os.path.exists(conversions_path) or not os.path.exists(index_path) or reparse or realign:
        logger.info('Parsing read, conversion, coverage information from BAM to' f'{conversions_path}, {coverage_path}')
        conversions_path, index_path, coverage_path = bam.parse_all_reads(
            STAR_result['bam'],
            conversions_path,
            index_path,
            coverage_path,
            n_threads=n_threads,
            temp_dir=temp_dir,
        )
    else:
        logger.info(
            'Skipping read and conversion parsing from BAM because files '
            'already exist. Use the `--reparse` flag to parse BAM alignments again.'
        )

    # Count conversions
    count_path = os.path.join(out_dir, constants.COUNT_FILENAME)
    if not os.path.exists(count_path) or recount or reparse or realign:
        logger.info(f'Counting conversions to {count_path}')
        count_path = conversions.count_conversions(
            conversions_path, index_path, count_path, quality=quality, n_threads=n_threads, temp_dir=temp_dir
        )
    else:
        logger.info(
            'Skipping conversion counting because files '
            'already exist. Use the `--recount` flag to count conversions again.'
        )
