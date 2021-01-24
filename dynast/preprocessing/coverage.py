import gzip
import logging
import os
import pickle
import shutil
import tempfile
from functools import partial

import pysam

from .. import utils

logger = logging.getLogger(__name__)


def calculate_coverage_contig(
    bam_path,
    counter,
    lock,
    contig,
    indices,
    umi_tag=None,
    barcode_tag=None,
    barcodes=None,
    temp_dir=None,
    update_every=50000,
    velocity=True
):
    """Calculate converage for a specific contig. This function is designed to
    be called as a separate process.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param contig: only reads that map to this contig will be processed
    :type contig: str
    :param indices: genomic positions to consider
    :type indices: list
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param update_every: update the counter every this many reads, defaults to `50000`
    :type update_every: int, optional
    :param velocity: whether or not velocities were assigned
    :type velocity: bool, optional

    :return: (coverage CSV path, index path)
    :rtype: (str, str)
    """
    coverage_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)

    # Index will be used for fast splitting, just like for conversions.csv
    index = []

    paired = {}
    coverage = {}

    required_tags = []
    if not velocity:
        required_tags.append('GX')
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    n = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam, \
        open(coverage_path, 'w') as coverage_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = coverage_out.tell()
            n += 1
            if n % update_every == 0:
                # Flush
                for key in sorted(list(coverage.keys())):
                    if key[1] < read.reference_start:
                        n_lines += 1
                        coverage_out.write(f'{key[0]},{contig},{key[1]},{coverage[key]}\n')
                        del coverage[key]

                lock.acquire()
                counter.value += update_every
                lock.release()

            if any(not read.has_tag(tag) for tag in required_tags):
                continue

            barcode = read.get_tag(barcode_tag) if barcode_tag else 'NA'
            if barcodes and barcode not in barcodes:
                continue

            read_id = read.query_name
            reference_positions = read.get_reference_positions()

            # For paired end reads, don't count overlap twice
            if read.is_paired:
                if read_id not in paired:
                    paired[read_id] = set(reference_positions)
                    continue

                reference_positions = list(paired[read_id].union(reference_positions))
                del paired[read_id]

            for genome_i in reference_positions:
                if genome_i in indices:
                    key = (barcode, genome_i)
                    coverage.setdefault(key, 0)
                    coverage[key] += 1

            if n_lines > 0:
                index.append((pos, n_lines))

        n_lines = 0
        pos = coverage_out.tell()
        for key in sorted(list(coverage.keys())):
            n_lines += 1
            coverage_out.write(f'{key[0]},{contig},{key[1]},{coverage[key]}\n')
        if n_lines > 0:
            index.append((pos, n_lines))
    lock.acquire()
    counter.value += n % update_every
    lock.release()
    index_path = utils.write_pickle(index, index_path, protocol=4)

    return coverage_path, index_path


def calculate_coverage(
    bam_path,
    conversions,
    coverage_path,
    index_path,
    umi_tag=None,
    barcode_tag=None,
    barcodes=None,
    n_threads=8,
    temp_dir=None,
    velocity=True,
):
    """Calculate coverage of each genomic position per barcode.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions: dictionary of contigs as keys and sets of genomic positions
                        as values that indicates positions where conversions were observed
    :type conversions: dictionary
    :param coverage_path: path to write coverage CSV
    :type coverage_path: str
    :param index_path: path to write index
    :type index_path: str
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param velocity: whether or not velocities were assigned
    :type velocity: bool, optional

    :return: (coverage CSV path, index path)
    :rtype: (str, str)
    """
    logger.debug(f'Extracting contigs from BAM {bam_path}')
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            calculate_coverage_contig,
            bam_path,
            counter,
            lock,
            umi_tag=umi_tag,
            barcode_tag=barcode_tag,
            barcodes=barcodes,
            temp_dir=tempfile.mkdtemp(dir=temp_dir),
            velocity=velocity
        ), [(contig, conversions.get(contig, set())) for contig in contigs]
    )
    pool.close()

    # Display progres bar
    logger.debug(f'Processing contigs {contigs} from BAM')
    utils.display_progress_with_counter(counter, n_reads, async_result)
    pool.join()

    # Combine csvs
    logger.debug(f'Writing coverage to {coverage_path}')
    pos = 0
    index = []
    with open(coverage_path, 'wb') as coverage_out:
        coverage_out.write(b'barcode,contig,genome_i,coverage\n')
        pos = coverage_out.tell()

        for coverage_part_path, index_part_path in async_result.get():
            with open(coverage_part_path, 'rb') as f:
                shutil.copyfileobj(f, coverage_out)

            # Parse index
            with gzip.open(index_part_path, 'rb') as f:
                index_part = pickle.load(f)
            for p, n in index_part:
                index.append((pos + p, n))
            pos += os.path.getsize(coverage_part_path)

    logger.debug(f'Writing coverage index to {index_path}')
    index_path = utils.write_pickle(index, index_path, protocol=4)

    return coverage_path, index_path
