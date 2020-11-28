import logging
import re
from functools import partial
from operator import truediv

import numpy as np
import pandas as pd

from .. import utils
from .index import split_index

logger = logging.getLogger(__name__)

CONVERSIONS_PARSER = re.compile(
    r'''^
    ([^,]*),
    (?P<barcode>[^,]*),
    ([^,]*?),?
    ([^,]*),
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    ([^,]*),
    ([^,]*),
    (?P<quality>[^,]*),
    ([^,]*),
    ([^,]*),
    ([^,]*),
    ([^,]*)\n
    $''', re.VERBOSE
)

COVERAGE_PARSER = re.compile(
    r'''^
    (?P<barcode>[^,]*),
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<coverage>[^,]*)\n
    $''', re.VERBOSE
)


def read_snps(snps_path, group_by=None):
    df = pd.read_csv(
        snps_path, dtype={
            'barcode': 'string',
            'contig': 'category',
            'genome_i': np.uint32,
        }
    )
    if group_by is None:
        return dict(df.groupby('contig').agg(set)['genome_i'])
    else:
        # TODO
        raise Exception()


def read_snp_csv(snp_csv):
    df = pd.read_csv(snp_csv, names=['contig', 'genome_i'])
    return dict(df.groupby('contig').agg(set)['genome_i'])


def extract_conversions_part(
    conversions_path, counter, lock, pos, n_lines, group_by=None, quality=27, update_every=10000
):
    conversions = {}
    with open(conversions_path, 'r') as f:
        f.seek(pos)

        for i in range(n_lines):
            line = f.readline()
            groupdict = CONVERSIONS_PARSER.match(line).groupdict()
            if int(groupdict['quality']) > quality:
                contig = groupdict['contig']
                genome_i = int(groupdict['genome_i'])

                if group_by is None:
                    conversions.setdefault(contig, {}).setdefault(genome_i, 0)
                    conversions[contig][genome_i] += 1
                else:
                    # TODO
                    pass
            if (i + 1) % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()
    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return conversions


def extract_conversions(conversions_path, index_path, group_by=None, quality=27, n_threads=8):
    logger.debug(f'Loading index {index_path} for {conversions_path}')
    index = utils.read_pickle(index_path)

    logger.debug(f'Splitting index into {n_threads} parts')
    parts = split_index(index, n=n_threads)

    logger.debug(f'Spawning {n_threads} processes')
    n_lines = sum(idx[1] for idx in index)
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            extract_conversions_part,
            conversions_path,
            counter,
            lock,
            group_by=group_by,
            quality=quality,
        ), parts
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(async_result, counter, n_lines)
    pool.join()

    logger.debug('Combining conversions')
    conversions = {}
    for conversions_part in async_result.get():
        conversions = utils.merge_dictionaries(conversions, conversions_part)

    return conversions


def extract_coverage_part(coverage_path, counter, lock, pos, n_lines, group_by=None, quality=27, update_every=10000):
    coverage = {}
    with open(coverage_path, 'r') as f:
        f.seek(pos)

        for i in range(n_lines):
            line = f.readline()
            groupdict = COVERAGE_PARSER.match(line).groupdict()

            contig = groupdict['contig']
            genome_i = int(groupdict['genome_i'])
            count = int(groupdict['coverage'])

            if group_by is None:
                coverage.setdefault(contig, {}).setdefault(genome_i, 0)
                coverage[contig][genome_i] += count
            else:
                # TODO
                pass
            if (i + 1) % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()
    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return coverage


def extract_coverage(coverage_path, index_path, group_by=None, quality=27, n_threads=8):
    logger.debug(f'Loading index {index_path} for {coverage_path}')
    index = utils.read_pickle(index_path)

    logger.debug(f'Splitting index into {n_threads} parts')
    parts = split_index(index, n=n_threads)

    logger.debug(f'Spawning {n_threads} processes')
    n_lines = sum(idx[1] for idx in index)
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            extract_coverage_part,
            coverage_path,
            counter,
            lock,
            group_by=group_by,
            quality=quality,
        ), parts
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(async_result, counter, n_lines)
    pool.join()

    logger.debug('Combining coverage')
    coverage = {}
    for coverage_part in async_result.get():
        coverage = utils.merge_dictionaries(coverage, coverage_part)

    return coverage


def detect_snps(
    conversions_path,
    conversions_index_path,
    coverage_path,
    coverage_index_path,
    snps_path,
    group_by=None,
    quality=27,
    threshold=0.5,
    n_threads=8,
):
    logger.info('Counting number of conversions for each genomic position')
    conversions = extract_conversions(
        conversions_path, conversions_index_path, group_by=group_by, quality=quality, n_threads=n_threads
    )

    logger.info('Counting coverage for each genomic position')
    coverage = extract_coverage(
        coverage_path, coverage_index_path, group_by=group_by, quality=quality, n_threads=n_threads
    )

    logger.info('Calculating fraction of conversions for each genomic position')
    fractions = utils.merge_dictionaries(conversions, coverage, f=truediv)

    logger.info(f'Writing detected SNPs to {snps_path}')
    with open(snps_path, 'w') as f:
        prefix = '' if group_by is None else f'{",".join(group_by)},'
        f.write(f'{prefix}contig,genome_i\n')
        for key, fraction in utils.flatten_dictionary(fractions):
            # If (# conversions) / (# coverage) is greater than a threshold,
            # consider this a SNP and write to CSV
            if fraction > threshold:
                f.write(f'{",".join(str(k) for k in key)}\n')

    return snps_path
