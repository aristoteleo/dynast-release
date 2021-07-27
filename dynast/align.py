import os
import tempfile

from . import config, constants, utils
from .logging import logger


def STAR_solo(
    fastqs,
    index_dir,
    out_dir,
    technology,
    whitelist_path=None,
    strand='forward',
    n_threads=8,
    temp_dir=None,
    nasc=False,
    overrides=None,
):
    """Align FASTQs with STARsolo.

    :param fastqs: list of path to FASTQs. Order matters -- STAR assumes the
                   UMI and barcode are in read 2
    :type fastqs: list
    :param index_dir: path to directory containing STAR index
    :type index_dir: str
    :param out_dir: path to directory to place STAR output
    :type out_dir: str
    :param technology: a `Technology` object defined in `technology.py`
    :type technology: collections.namedtuple
    :param whitelist_path: path to textfile containing barcode whitelist,
                           defaults to `None`
    :type whitelist_path: str, optional
    :param strand: strandedness of the sequencing protocol, defaults to `forward`,
                   may be one of the following: `forward`, `reverse`, `unstranded`
    :type strand: str, optional
    :param n_threads: number of threads to use, defaults to `8`
    :type n_threads: int, optional
    :param temp_dir: STAR temporary directory, defaults to `None`, which
                     uses the system temporary directory
    :type temp_dir: str, optional
    :param nasc: whether or not to use STAR configuration used in NASC-seq pipeline,
                 defaults to `False`
    :type nasc: bool, optional
    :param overrides: STAR command-line argument overrides, defaults to `None`
    :type overrides: dictionary, optional

    :return: dictionary containing output files
    :rtype: dict
    """
    logger.info('Aligning the following FASTQs with STAR')
    for fastq in fastqs:
        logger.info((' ' * 8) + fastq)

    # out_dir must end with a directory separator
    if not out_dir.endswith(('/', '\\')):
        out_dir += os.path.sep

    command = [utils.get_STAR_binary_path()]
    arguments = config.STAR_ARGUMENTS
    arguments = utils.combine_arguments(arguments, technology.chemistry.to_starsolo_arguments())
    if technology.additional_args:
        arguments = utils.combine_arguments(arguments, technology.additional_args)
    if technology.name != 'smartseq':
        arguments = utils.combine_arguments(arguments, config.STAR_SOLO_ARGUMENTS)
        if whitelist_path:
            arguments = utils.combine_arguments(arguments, {'--soloCBwhitelist': whitelist_path})

        # Input FASTQs must be plaintext
        plaintext_fastqs = []
        for fastq in fastqs:
            if fastq.endswith('.gz'):
                plaintext_path = utils.mkstemp(dir=temp_dir)
                logger.warning(f'Decompressing {fastq} to {plaintext_path} because STAR requires plaintext FASTQs')
                utils.decompress_gzip(fastq, plaintext_path)
            else:
                plaintext_path = fastq
            plaintext_fastqs.append(plaintext_path)
        arguments['--readFilesIn'] = plaintext_fastqs
    else:
        # STAR requires FIFO file support when running in smartseq mode
        fifo_path = utils.mkstemp(dir=temp_dir, delete=True)
        try:
            os.mkfifo(fifo_path)
        except OSError:
            raise utils.UnsupportedOSException(
                f'The filesystem at {temp_dir} does not support FIFO files. '
                'STAR uses FIFO files to run alignment of Smart-seq files.'
            )

        if nasc:
            arguments = utils.combine_arguments(arguments, config.NASC_ARGUMENTS)

        manifest_path = utils.mkstemp(dir=temp_dir)
        with open(fastqs[0], 'r') as f, open(manifest_path, 'w') as out:
            for line in f:
                if line.isspace():
                    continue

                cell_id, fastq_1, fastq_2 = line.strip().split(',')
                if fastq_1.endswith('.gz'):
                    plaintext_path = utils.mkstemp(dir=temp_dir)
                    logger.warning(
                        f'Decompressing {fastq_1} to {plaintext_path} because STAR requires plaintext FASTQs'
                    )
                    utils.decompress_gzip(fastq_1, plaintext_path)
                    fastq_1 = plaintext_path
                if fastq_2.endswith('.gz'):
                    plaintext_path = utils.mkstemp(dir=temp_dir)
                    logger.warning(
                        f'Decompressing {fastq_2} to {plaintext_path} because STAR requires plaintext FASTQs'
                    )
                    utils.decompress_gzip(fastq_2, plaintext_path)
                    fastq_2 = plaintext_path

                out.write(f'{fastq_1}\t{fastq_2}\t{cell_id}\n')
        arguments['--readFilesManifest'] = manifest_path

    # Attempt to increase NOFILE limit to maximum possible, and use this to
    # set n_bins to the maximum possible, as higher n_bin requires less
    # memory to sort the BAM.
    current = utils.get_file_descriptor_limit()
    maximum = min(utils.get_max_file_descriptor_limit(), 100000)
    arguments = utils.combine_arguments(
        arguments, {
            '--soloStrand': strand.capitalize(),
            '--genomeDir': index_dir,
            '--runThreadN': n_threads,
            '--outFileNamePrefix': out_dir,
            '--outTmpDir': os.path.join(temp_dir, f'{tempfile.gettempprefix()}{next(tempfile._get_candidate_names())}'),
            '--outBAMsortingBinsN': max((maximum // n_threads) - 10, 50)
        }
    )
    arguments.update(overrides or {})

    logger.debug(f'Increasing maximum number of open file descriptors from {current} to {maximum}')
    with utils.increase_file_descriptor_limit(maximum):
        logger.info('Starting alignment')
        command += utils.arguments_to_list(arguments)
        utils.run_executable(command)

    solo_dir = os.path.join(out_dir, constants.STAR_SOLO_DIR)
    gene_dir = os.path.join(solo_dir, constants.STAR_GENE_DIR)
    raw_gene_dir = os.path.join(gene_dir, constants.STAR_RAW_DIR)
    filtered_gene_dir = os.path.join(gene_dir, constants.STAR_FILTERED_DIR)

    result = {
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
    }
    return result


@logger.namespaced('align')
def align(
    fastqs,
    index_dir,
    out_dir,
    technology,
    whitelist_path=None,
    strand='forward',
    n_threads=8,
    temp_dir=None,
    nasc=False,
    overrides=None,
):
    os.makedirs(out_dir, exist_ok=True)

    # Check memory.
    available_memory = utils.get_available_memory()
    if available_memory < config.RECOMMENDED_MEMORY:
        logger.warning(
            f'There is only {available_memory / (1024 ** 3):.2f} GB of free memory on the machine. '
            f'It is highly recommended to have at least {config.RECOMMENDED_MEMORY // (1024 ** 3)} GB '
            'free when running dynast. Continuing may cause dynast to crash with an out-of-memory error.'
        )

    # If whitelist was not provided but one is available, decompress into output
    # directory.
    if whitelist_path is None and technology.chemistry.has_whitelist:
        whitelist_path = os.path.join(out_dir, f'{technology.name}_whitelist.txt')
        logger.info(f'Copying prepackaged whitelist for technology {technology.name} to {whitelist_path}')
        utils.decompress_gzip(technology.chemistry.whitelist_path, whitelist_path)

    STAR_solo_dir = os.path.join(out_dir, constants.STAR_SOLO_DIR)
    STAR_gene_dir = os.path.join(STAR_solo_dir, constants.STAR_GENE_DIR)
    STAR_raw_gene_dir = os.path.join(STAR_gene_dir, constants.STAR_RAW_DIR)
    STAR_filtered_gene_dir = os.path.join(STAR_gene_dir, constants.STAR_FILTERED_DIR)
    STAR_velocyto_dir = os.path.join(STAR_solo_dir, constants.STAR_VELOCYTO_DIR, constants.STAR_RAW_DIR)
    # Check if these files exist. If they do, we can skip alignment.
    STAR_result = {
        'bam': os.path.join(out_dir, constants.STAR_BAM_FILENAME),
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
    skip = utils.all_exists(*STAR_required)
    if not skip:
        logger.info(f'STAR binary found at {utils.get_STAR_binary_path()}')
        STAR_result = STAR_solo(
            fastqs,
            index_dir,
            out_dir,
            technology,
            whitelist_path=whitelist_path,
            strand=strand,
            n_threads=n_threads,
            temp_dir=temp_dir,
            nasc=nasc,
            overrides=overrides,
        )
    else:
        logger.info('Alignment files already exist. Provide `--overwrite` to overwrite.')
