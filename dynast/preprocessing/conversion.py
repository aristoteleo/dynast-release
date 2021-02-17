import re
import shutil
import tempfile
from functools import partial

import numpy as np
import pandas as pd
from scipy import sparse

from .. import utils
from ..logging import logger

CONVERSIONS_PARSER = re.compile(
    r'''^
    (?P<read_id>[^,]*),
    (?P<barcode>[^,]*),
    (?P<umi>[^,]*),
    (?P<GX>[^,]*),
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<original>[^,]*),
    (?P<converted>[^,]*),
    (?P<quality>[^,]*),
    (?P<A>[^,]*),
    (?P<C>[^,]*),
    (?P<G>[^,]*),
    (?P<T>[^,]*),
    (?P<velocity>[^,]*),
    (?P<transcriptome>[^,]*)\n
    $''', re.VERBOSE
)

NO_CONVERSIONS_PARSER = re.compile(
    r'''^
    (?P<read_id>[^,]*),
    (?P<barcode>[^,]*),
    (?P<umi>[^,]*),
    (?P<GX>[^,]*),
    (?P<A>[^,]*),
    (?P<C>[^,]*),
    (?P<G>[^,]*),
    (?P<T>[^,]*),
    (?P<velocity>[^,]*),
    (?P<transcriptome>[^,]*)\n
    $''', re.VERBOSE
)

CONVERSION_IDX = {
    ('A', 'C'): 0,
    ('A', 'G'): 1,
    ('A', 'T'): 2,
    ('C', 'A'): 3,
    ('C', 'G'): 4,
    ('C', 'T'): 5,
    ('G', 'A'): 6,
    ('G', 'C'): 7,
    ('G', 'T'): 8,
    ('T', 'A'): 9,
    ('T', 'C'): 10,
    ('T', 'G'): 11,
}
BASE_IDX = {
    'A': 12,
    'C': 13,
    'G': 14,
    'T': 15,
}
CONVERSION_COLUMNS = [''.join(pair) for pair in sorted(CONVERSION_IDX.keys())]
BASE_COLUMNS = sorted(BASE_IDX.keys())
COLUMNS = CONVERSION_COLUMNS + BASE_COLUMNS


def read_counts(counts_path, *args, **kwargs):
    """Read counts CSV as a pandas dataframe.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    :param counts_path: path to CSV
    :type counts_path: str

    :return: counts dataframe
    :rtype: pandas.DataFrame
    """
    dtypes = {
        'read_id': 'string',
        'barcode': 'string',
        'umi': 'string',
        'GX': 'string',
        'velocity': 'category',
        'transcriptome': bool,
        **{column: np.uint8
           for column in COLUMNS}
    }
    return pd.read_csv(counts_path, dtype=dtypes, na_filter=False, *args, **kwargs)


def complement_counts(df_counts, gene_infos):
    """Complement the counts in the counts dataframe according to gene strand.

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame
    :param gene_infos: dictionary containing gene information, as returned by
                       `preprocessing.gtf.parse_gtf`
    :type gene_infos: dictionary

    :return: counts dataframe with counts complemented for reads mapping to genes on the reverse strand
    :rtype: pandas.DataFrame
    """
    df_counts['strand'] = df_counts['GX'].map(lambda gx: gene_infos[gx]['strand'])

    columns = ['barcode', 'GX', 'velocity', 'transcriptome'] + COLUMNS
    df_forward = df_counts[df_counts.strand == '+'][columns]
    df_reverse = df_counts[df_counts.strand == '-'][columns]

    df_reverse.columns = ['barcode', 'GX', 'velocity', 'transcriptome'] + CONVERSION_COLUMNS[::-1] + BASE_COLUMNS[::-1]
    df_reverse = df_reverse[columns]

    return pd.concat((df_forward, df_reverse)).reset_index()


def drop_multimappers(df_counts):
    """Drop multimappings that have the same read ID where
    * some map to the transcriptome while some do not -- drop non-transcriptome alignments
    * none map to the transcriptome AND aligned to multiple genes -- drop all
    * none map to the transcriptome AND assigned multiple velocity types -- set to ambiguous

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame

    :return: counts dataframe with multimappers appropriately filtered
    :rtype: pandas.DataFrame
    """
    # Multimapping reads
    duplicated_mask = df_counts.duplicated('read_id', keep=False)
    if not any(duplicated_mask):
        return df_counts

    df_multi = df_counts[duplicated_mask]

    filtered = []
    for read_id, df_read in df_multi.groupby('read_id'):
        transcriptome = list(df_read['transcriptome'])
        # Rule 1
        if True in transcriptome and False in transcriptome:
            filtered.append(df_read[transcriptome])
        # Rule 2, 3
        elif all(~df_read['transcriptome']):
            # Rule 2
            if len(df_read['GX'].unique()) > 1:
                continue
            # Rule 3
            elif len(df_read['velocity'].unique()) > 1:
                df_read = df_read.copy()
                df_read['velocity'] = 'ambiguous'
                filtered.append(df_read.drop_duplicates('read_id', keep='first'))
            else:
                filtered.append(df_read.drop_duplicates('read_id', keep='first'))

    return df_counts[~duplicated_mask].append(
        pd.concat(filtered, ignore_index=True), ignore_index=True
    ) if filtered else df_counts[~duplicated_mask]


def deduplicate_counts(df_counts, conversions=None):
    """Deduplicate counts based on barcode, UMI, and gene.

    The order of priority is the following.
    1. Reads that align to the transcriptome (exon only)
    2. If `conversions` is provided, reads that have a larger sum of such conversions
       If `conversions` is not provided, reads that have larger sum of all conversions

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame
    :param conversions: conversions to prioritize, defaults to `None`
    :type conversions: list, optional

    :return: deduplicated counts dataframe
    :rtype: pandas.DataFrame
    """
    df_counts['conversion_sum'] = df_counts[conversions or CONVERSION_COLUMNS].sum(axis=1)

    # Sort by transcriptome last, longest alignment last, least conversion first
    df_sorted = df_counts.sort_values(['transcriptome', 'conversion_sum']).drop(columns='conversion_sum')

    # Always select transcriptome read if there are duplicates
    df_deduplicated = df_sorted[~df_sorted.duplicated(subset=['barcode', 'umi', 'GX'], keep='last')].sort_values(
        'barcode'
    ).reset_index(drop=True)
    return df_deduplicated


def split_counts_by_velocity(df_counts):
    """Split the given counts dataframe by the `velocity` column.

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame

    :return: dictionary containing `velocity` column values as keys and the
             subset dataframe as values
    :rtype: dictionary
    """
    dfs = {}
    for velocity, df_part in df_counts.groupby('velocity'):
        dfs[velocity] = df_part.reset_index(drop=True)
    logger.debug(f'Found the following velocity assignments: {", ".join(dfs.keys())}')
    return dfs


def counts_to_matrix(df_counts, barcodes, features, barcode_column='barcode', feature_column='GX'):
    """Convert a counts dataframe to a sparse counts matrix.

    Counts are assumed to be appropriately deduplicated.

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame
    :param barcodes: list of barcodes that will map to the rows
    :type barcodes: list
    :param features: list of features (i.e. genes) that will map to the columns
    :type features: list
    :param barcode_column: column in counts dataframe to use as barcodes, defaults to `barcode`
    :type barcode_column: str
    :param feature_column: column in counts dataframe to use as features, defaults to `GX`
    :type feature_column: str

    :return: sparse counts matrix
    :rtype: scipy.sparse.csrmatrix
    """
    # Transform to index for fast lookup
    barcode_indices = {barcode: i for i, barcode in enumerate(barcodes)}
    feature_indices = {feature: i for i, feature in enumerate(features)}

    matrix = sparse.lil_matrix((len(barcodes), len(features)), dtype=np.uint32)
    for (barcode, feature), count in df_counts.groupby([barcode_column, feature_column]).size().items():
        matrix[barcode_indices[barcode], feature_indices[feature]] = count

    return matrix.tocsr()


def split_counts(df_counts, barcodes, features, barcode_column='barcode', feature_column='GX', conversions=['TC']):
    """Split counts dataframe into two count matrices by a column.

    :param df_counts: counts dataframe
    :type df_counts: pandas.DataFrame
    :param barcodes: list of barcodes that will map to the rows
    :type barcodes: list
    :param features: list of features (i.e. genes) that will map to the columns
    :type features: list
    :param barcode_column: column in counts dataframe to use as barcodes, defaults to `barcode`
    :type barcode_column: str, optional
    :param feature_column: column in counts dataframe to use as features, defaults to `GX`
    :type feature_column: str, optional
    :param conversions: conversion(s) in question, defaults to `['TC']`
    :type conversions: list, optional

    :return: (count matrix of `conversion==0`, count matrix of `conversion>0`)
    :rtype: (scipy.sparse.csrmatrix, scipy.sparse.csrmatrix)
    """
    matrix_unlabeled = counts_to_matrix(
        df_counts[(df_counts[conversions] == 0).all(axis=1)],
        barcodes,
        features,
        barcode_column=barcode_column,
        feature_column=feature_column
    )
    matrix_labeled = counts_to_matrix(
        df_counts[(df_counts[conversions] > 0).any(axis=1)],
        barcodes,
        features,
        barcode_column=barcode_column,
        feature_column=feature_column
    )
    return matrix_unlabeled, matrix_labeled


def count_no_conversions_part(
    no_conversions_path,
    counter,
    lock,
    pos,
    n_lines,
    barcodes=None,
    temp_dir=None,
    update_every=10000,
):
    """Count reads that have no conversion.

    :param no_conversions_path: no conversions CSV path
    :type no_conversions_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param pos: file handle position at which to start reading the conversions CSV
    :type pos: int
    :param n_lines: number of lines to parse from the conversions CSV, starting
                    from position `pos`
    :type n_lines: int
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `5000`
    :type update_every: int, optional

    :return: path to temporary counts CSV
    :rtype: str
    """
    count_path = utils.mkstemp(dir=temp_dir)
    with open(no_conversions_path, 'r') as f, open(count_path, 'w') as out:
        f.seek(pos)
        for i in range(n_lines):
            if (i + 1) % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            line = f.readline()
            groups = NO_CONVERSIONS_PARSER.match(line).groupdict()
            if barcodes and groups['barcode'] not in barcodes:
                continue
            out.write(
                f'{groups["read_id"]},{groups["barcode"]},{groups["umi"]},{groups["GX"]},'
                f'{",".join(groups.get(key, "0") for key in COLUMNS)},{groups["velocity"]},{groups["transcriptome"]}\n'
            )
    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return count_path


def count_conversions_part(
    conversions_path,
    counter,
    lock,
    pos,
    n_lines,
    barcodes=None,
    snps=None,
    quality=27,
    temp_dir=None,
    update_every=10000,
):
    """Count the number of conversions of each read per barcode and gene, along with
    the total nucleotide content of the region each read mapped to, also per barcode
    and gene. This function is used exclusively for multiprocessing.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param pos: file handle position at which to start reading the conversions CSV
    :type pos: int
    :param n_lines: number of lines to parse from the conversions CSV, starting
                    from position `pos`
    :type n_lines: int
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param snps: dictionary of contig as keys and list of genomic positions as
                 values that indicate SNP locations, defaults to `None`
    :type snps: dictionary, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `10000`
    :type update_every: int, optional

    :return: path to temporary counts CSV
    :rtype: tuple
    """

    def is_snp(g):
        if not snps:
            return False
        return int(g['genome_i']) in snps.get(g['contig'], set())

    count_path = utils.mkstemp(dir=temp_dir)

    counts = None
    read_id = None
    count_base = True
    with open(conversions_path, 'r') as f, open(count_path, 'w') as out:
        f.seek(pos)

        groups = None
        prev_groups = None
        for i in range(n_lines):
            if (i + 1) % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            line = f.readline()
            prev_groups = groups
            groups = CONVERSIONS_PARSER.match(line).groupdict()

            if read_id != groups['read_id']:
                if read_id is not None and (not barcodes or prev_groups['barcode'] in barcodes):
                    out.write(
                        f'{prev_groups["read_id"]},{prev_groups["barcode"]},{prev_groups["umi"]},{prev_groups["GX"]},'
                        f'{",".join(str(c) for c in counts)},{prev_groups["velocity"]},{prev_groups["transcriptome"]}\n'
                    )
                counts = [0] * (len(CONVERSION_IDX) + len(BASE_IDX))
                read_id = groups['read_id']
                count_base = True
            if int(groups['quality']) > quality and not is_snp(groups):
                counts[CONVERSION_IDX[(groups['original'], groups['converted'])]] += 1
            if count_base:
                for base, j in BASE_IDX.items():
                    counts[j] = int(groups[base])
                count_base = False

        # Add last record
        if read_id is not None and (not barcodes or groups['barcode'] in barcodes):
            out.write(
                f'{groups["read_id"]},{groups["barcode"]},{groups["umi"]},{groups["GX"]},'
                f'{",".join(str(c) for c in counts)},{groups["velocity"]},{groups["transcriptome"]}\n'
            )

    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return count_path


def count_conversions(
    conversions_path,
    index_path,
    no_conversions_path,
    no_index_path,
    counts_path,
    barcodes=None,
    snps=None,
    quality=27,
    conversions=None,
    n_threads=8,
    temp_dir=None
):
    """Count the number of conversions of each read per barcode and gene, along with
    the total nucleotide content of the region each read mapped to, also per barcode.
    When a duplicate UMI for a barcode is observed, the read with the greatest
    number of conversions is selected.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param no_conversions_path: path to output information about reads that do not have any conversions
    :type no_conversions_path: str
    :param no_index_path: path to no conversions index
    :type no_index_path: str
    :param counts_path: path to write counts CSV
    :param counts_path: str
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param snps: dictionary of contig as keys and list of genomic positions as
                 values that indicate SNP locations, defaults to `None`
    :type snps: dictionary, optional
    :param conversions: conversions to prioritize when deduplicating only applicable
                        for UMI technologies, defaults to `None`
    :type conversions: list, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional

    :return: path to counts CSV
    :rtype: str
    """
    # Load index
    logger.debug(f'Loading index {index_path} for {conversions_path}')
    idx = utils.read_pickle(index_path)
    no_idx = utils.read_pickle(no_index_path)

    # Split index into n contiguous pieces
    logger.debug(f'Splitting index into {n_threads} parts')
    parts = utils.split_index(idx, n=n_threads)
    no_parts = []
    for i in range(0, len(no_idx), (len(no_idx) // n_threads) + 1):
        no_parts.append((no_idx[i], min((len(no_idx) // n_threads) + 1, len(no_idx[i:]))))

    # Parse each part in a different process
    logger.debug(f'Spawning {n_threads} processes')
    n_lines = sum(i[1] for i in idx) + len(no_idx)
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            count_conversions_part,
            conversions_path,
            counter,
            lock,
            barcodes=barcodes,
            snps=snps,
            quality=quality,
            temp_dir=tempfile.mkdtemp(dir=temp_dir)
        ), parts
    )
    no_async_result = pool.starmap_async(
        partial(
            count_no_conversions_part,
            no_conversions_path,
            counter,
            lock,
            barcodes=barcodes,
            temp_dir=tempfile.mkdtemp(dir=temp_dir)
        ), no_parts
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, n_lines, async_result, no_async_result)
    pool.join()

    # Combine csvs
    combined_path = utils.mkstemp(dir=temp_dir)
    logger.debug(f'Combining intermediate parts to {combined_path}')
    with open(combined_path, 'wb') as out:
        out.write(f'read_id,barcode,umi,GX,{",".join(COLUMNS)},velocity,transcriptome\n'.encode())
        for counts_part_path in async_result.get():
            with open(counts_part_path, 'rb') as f:
                shutil.copyfileobj(f, out)
        for counts_part_path in no_async_result.get():
            with open(counts_part_path, 'rb') as f:
                shutil.copyfileobj(f, out)

    # Deduplicate using UMIs if possible, otherwise use read id
    df_counts = read_counts(combined_path)
    if all(df_counts['umi'] != 'NA'):
        logger.debug(f'Deduplicating reads based on barcode and UMI to {counts_path}')
        deduplicate_counts(df_counts, conversions=conversions).drop(columns='read_id').to_csv(counts_path, index=False)
    else:
        logger.debug(f'Filtering multimappers to {counts_path}')
        drop_multimappers(df_counts).drop(columns='read_id').to_csv(counts_path, index=False)

    return counts_path
