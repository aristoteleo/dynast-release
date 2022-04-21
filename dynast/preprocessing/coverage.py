import re
from collections import Counter
from functools import partial

import ngs_tools as ngs
import pysam

from . import bam
from .. import utils
from ..logging import logger

COVERAGE_PARSER = re.compile(
    r'''^
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<coverage>[^,]*)\n
    $''', re.VERBOSE
)


def read_coverage(coverage_path):
    """Read coverage CSV as a dictionary.

    :param coverage_path: path to coverage CSV
    :type coverage_path: str

    :return: coverage as a nested dictionary
    :rtype: dict
    """
    coverage = {}
    with open(coverage_path, 'r') as f:
        for i, line in enumerate(f):
            groups = COVERAGE_PARSER.match(line).groupdict()
            try:
                contig, genome_i, count = groups['contig'], int(groups['genome_i']), int(groups['coverage'])
            except ValueError as e:
                if i == 0:
                    continue
                raise e
            coverage.setdefault(contig, {})[genome_i] = count
    return coverage


def calculate_coverage_contig(
    counter,
    lock,
    bam_path,
    contig,
    indices,
    alignments=None,
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    barcodes=None,
    temp_dir=None,
    update_every=10000,
    velocity=True
):
    """Calculate converage for a specific contig. This function is designed to
    be called as a separate process.

    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param contig: only reads that map to this contig will be processed
    :type contig: str
    :param indices: genomic positions to consider
    :type indices: list
    :param alignments: set of (read_id, alignment_index) tuples to process. All
        alignments are processed if this option is not provided.
    :type alignments: set, optional
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
    :type gene_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `10000`
    :type update_every: int, optional
    :param velocity: whether or not velocities were assigned
    :type velocity: bool, optional

    :return: coverag
    :rtype: dict
    """
    coverage_path = utils.mkstemp(dir=temp_dir)
    paired = {}
    coverage = Counter()

    required_tags = []
    if not velocity:
        required_tags.append(gene_tag)
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(coverage_path, 'w') as coverage_out:
        for read in bam.fetch(contig):
            n += 1
            if n % update_every == 0:
                # NOTE: dictionary keys are sorted by insertion order
                for genome_i in list(coverage.keys()):
                    if genome_i < read.reference_start:
                        coverage_out.write(f'{contig},{genome_i},{coverage[genome_i]}\n')
                        del coverage[genome_i]
                    else:
                        break
                lock.acquire()
                counter.value += update_every
                lock.release()

            if any(not read.has_tag(tag) for tag in required_tags) or read.is_duplicate or read.is_unmapped:
                continue

            barcode = read.get_tag(barcode_tag) if barcode_tag else 'NA'
            if barcodes and barcode not in barcodes:
                continue

            read_id = read.query_name
            alignment_index = read.get_tag('HI')
            key = (read_id, alignment_index)
            if alignments and key not in alignments:
                continue

            reference_positions = read.get_reference_positions()

            # For paired end reads, don't count overlap twice
            if read.is_paired:
                if key not in paired:
                    paired[key] = reference_positions
                    continue

                reference_positions = list(set(reference_positions + paired[key]))
                del paired[key]

            for genome_i in reference_positions:
                if genome_i in indices:
                    coverage[genome_i] += 1

        for genome_i, cover in coverage.items():
            coverage_out.write(f'{contig},{genome_i},{cover}\n')

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    return coverage_path


def calculate_coverage(
    bam_path,
    conversions,
    coverage_path,
    alignments=None,
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    barcodes=None,
    n_threads=8,
    temp_dir=None,
    velocity=True,
    splits=None,
):
    """Calculate coverage of each genomic position per barcode.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions: dictionary of contigs as keys and sets of genomic positions
                        as values that indicates positions where conversions were observed
    :type conversions: dictionary
    :param coverage_path: path to write coverage CSV
    :type coverage_path: str
    :param alignments: set of (read_id, alignment_index) tuples to process. All
        alignments are processed if this option is not provided.
    :type alignments: set, optional
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
    :type gene_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param velocity: whether or not velocities were assigned
    :type velocity: bool, optional
    :param splits: list containing tuples of (pre-split BAM paths, number of reads),
                   defaults to `None`
    :type splits: list, optional

    :return: coverage CSV path
    :rtype: str
    """
    if splits is None:
        splits = bam.split_bam(
            bam_path,
            min(n_threads * 8,
                utils.get_file_descriptor_limit() - 10),
            n_threads=n_threads,
            temp_dir=temp_dir
        ) if n_threads > 1 else [(bam_path, ngs.bam.count_bam(bam_path))]

    # Add argument per contig
    args = []
    for path, _ in splits:
        with pysam.AlignmentFile(path, 'rb') as f:
            for index in f.get_index_statistics():
                if index.total > 0:
                    args.append((path, index.contig, conversions.get(index.contig, set())))

    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            calculate_coverage_contig,
            counter,
            lock,
            alignments=alignments,
            umi_tag=umi_tag,
            barcode_tag=barcode_tag,
            gene_tag=gene_tag,
            barcodes=barcodes,
            velocity=velocity
        ), args
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, sum(value[1] for value in splits), async_result)
    pool.join()

    logger.debug('Combining results')
    coverage = {}
    for coverage_part_path in async_result.get():
        cov = read_coverage(coverage_part_path)
        coverage = utils.merge_dictionaries(coverage, cov)

    logger.debug(f'Writing coverage to {coverage_path}')
    with open(coverage_path, 'w') as coverage_out:
        coverage_out.write('contig,genome_i,coverage\n')
        for (contig, genome_i), cover in utils.flatten_dictionary(coverage):
            coverage_out.write(f'{contig},{genome_i},{cover}\n')

    return coverage_path
