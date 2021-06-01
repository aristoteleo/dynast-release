import re
from functools import partial
from operator import truediv

import numpy as np
import pandas as pd

from .. import utils
from ..logging import logger
from .conversion import CONVERSIONS_PARSER

COVERAGE_PARSER = re.compile(
    r'''^
    (?P<barcode>[^,]*),
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<coverage>[^,]*)\n
    $''', re.VERBOSE
)


def read_snps(snps_path):
    """Read SNPs CSV as a dictionary

    :param snps_path: path to SNPs CSV
    :type snps_path: str

    :return: dictionary of contigs as keys and sets of genomic positions with SNPs as values
    :rtype: dictionary
    """
    df = pd.read_csv(
        snps_path, dtype={
            'barcode': 'category',
            'contig': 'category',
            'genome_i': np.uint32,
        }
    )
    return dict(df.groupby('contig', sort=False, observed=True).agg(set)['genome_i'])


def read_snp_csv(snp_csv):
    """Read a user-provided SNPs CSV

    :param snp_csv: path to SNPs CSV
    :type snp_csv: str

    :return: dictionary of contigs as keys and sets of genomic positions with SNPs as values
    :rtype: dictionary
    """
    # Check if header exists
    header = False
    with open(snp_csv, 'r') as f:
        line = f.readline().strip()
        try:
            int(line.split(',')[-1])
        except ValueError:
            header = True

    df = pd.read_csv(
        snp_csv,
        names=['contig', 'genome_i'],
        skiprows=1 if header else None,
        dtype={
            'contig': 'category',
            'genome_i': np.uint32
        }
    )
    return dict(df.groupby('contig', sort=False, observed=True).agg(set)['genome_i'])


def extract_conversions_part(conversions_path, counter, lock, index, alignments=None, quality=27, update_every=10000):
    """Extract number of conversions for every genomic position.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param index: list of (file position, number of lines) tuples to process
    :type index: list
    :param alignments: set of (read_id, alignment_index) tuples to process. All
        alignments are processed if this option is not provided.
    :type alignments: set, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param update_every: update the counter every this many reads, defaults to `10000`
    :type update_every: int, optional

    :return: nested dictionary that contains number of conversions for each contig and position
    :rtype: dictionary
    """
    conversions = {}
    n = 0
    with open(conversions_path, 'r') as f:
        for pos, n_lines, pos2 in index:
            f.seek(pos)
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            for _ in range(n_lines):
                line = f.readline()
                groups = CONVERSIONS_PARSER.match(line).groupdict()
                key = (groups['read_id'], int(groups['index']))
                if alignments and key not in alignments:
                    break

                if int(groups['quality']) > quality:
                    contig = groups['contig']
                    genome_i = int(groups['genome_i'])

                    conversions.setdefault(contig, {}).setdefault(genome_i, 0)
                    conversions[contig][genome_i] += 1
    lock.acquire()
    counter.value += n % update_every
    lock.release()

    return conversions


def extract_conversions(conversions_path, index_path, alignments=None, quality=27, n_threads=8):
    """Wrapper around `extract_conversions_part` that works in parallel.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param alignments: set of (read_id, alignment_index) tuples to process. All
        alignments are processed if this option is not provided.
    :type alignments: set, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: nested dictionary that contains number of conversions for each contig and position
    :rtype: dictionary
    """
    logger.debug(f'Loading index {index_path} for {conversions_path}')
    index = utils.read_pickle(index_path)

    logger.debug(f'Splitting index into {n_threads} parts')
    parts = utils.split_index(index, n=n_threads)

    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            extract_conversions_part,
            conversions_path,
            counter,
            lock,
            alignments=alignments,
            quality=quality,
        ), [(part,) for part in parts]
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, len(index), async_result)
    pool.join()

    logger.debug('Combining conversions')
    conversions = {}
    for conversions_part in async_result.get():
        conversions = utils.merge_dictionaries(conversions, conversions_part)

    return conversions


def extract_coverage_part(coverage_path, counter, lock, index, update_every=5000):
    """Extract coverage for every genomic position.

    :param coverage_path: path to coverage CSV
    :type coverage_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param index: list of (file position, number of lines) tuples to process
    :type index: list
    :param update_every: update the counter every this many reads, defaults to `5000`
    :type update_every: int, optional

    :return: nested dictionary that contains number of conversions for each contig and position
    :rtype: dictionary
    """
    coverage = {}
    n = 0
    with open(coverage_path, 'r') as f:
        for tup in index:
            pos = tup[0]
            n_lines = tup[1]
            f.seek(pos)
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            for i in range(n_lines):
                line = f.readline()
                groupdict = COVERAGE_PARSER.match(line).groupdict()

                contig = groupdict['contig']
                genome_i = int(groupdict['genome_i'])
                count = int(groupdict['coverage'])

                coverage.setdefault(contig, {}).setdefault(genome_i, 0)
                coverage[contig][genome_i] += count

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    return coverage


def extract_coverage(coverage_path, index_path, n_threads=8):
    """Wrapper around `extract_coverage_part` that works in parallel.

    :param coverage_path: path to coverage CSV
    :type coverage_path: str
    :param index_path: path to coverage index
    :type index_path: str
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: nested dictionary that contains number of conversions for each contig and position
    :rtype: dictionary
    """
    logger.debug(f'Loading index {index_path} for {coverage_path}')
    index = utils.read_pickle(index_path)

    logger.debug(f'Splitting index into {n_threads} parts')
    parts = utils.split_index(index, n=n_threads)

    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            extract_coverage_part,
            coverage_path,
            counter,
            lock,
        ), [(part,) for part in parts]
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, len(index), async_result)
    pool.join()

    logger.debug('Combining coverage')
    coverage = {}
    for coverage_part in async_result.get():
        coverage = utils.merge_dictionaries(coverage, coverage_part)

    return coverage


def detect_snps(
    conversions_path,
    index_path,
    coverage_path,
    coverage_index_path,
    snps_path,
    alignments=None,
    quality=27,
    threshold=0.5,
    n_threads=8,
):
    """Detect SNPs.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param coverage_path: path to coverage CSV
    :type coverage_path: str
    :param coverage_index_path: path to coverage index
    :type coverage_index_path: str
    :param snps_path: path to output SNPs
    :type snps_path: str
    :param alignments: set of (read_id, alignment_index) tuples to process. All
        alignments are processed if this option is not provided.
    :type alignments: set, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param threshold: positions with conversions / coverage > threshold will be
                      considered as SNPs, defaults to `0.5`
    :type threshold: float, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    """
    logger.debug('Counting number of conversions for each genomic position')
    conversions = extract_conversions(
        conversions_path, index_path, alignments=alignments, quality=quality, n_threads=n_threads
    )

    logger.debug('Counting coverage for each genomic position')
    coverage = extract_coverage(coverage_path, coverage_index_path, n_threads=n_threads)

    logger.debug('Calculating fraction of conversions for each genomic position')
    fractions = utils.merge_dictionaries(conversions, coverage, f=truediv)

    logger.debug(f'Writing detected SNPs to {snps_path}')
    with open(snps_path, 'w') as f:
        f.write('contig,genome_i\n')
        for (contig, genome_i), fraction in utils.flatten_dictionary(fractions):
            # If (# conversions) / (# coverage) is greater than a threshold,
            # consider this a SNP and write to CSV
            if fraction > threshold:
                f.write(f'{contig},{genome_i}\n')

    return snps_path
