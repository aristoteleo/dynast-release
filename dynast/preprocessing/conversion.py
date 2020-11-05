import logging
import re
import shutil
import tempfile
from functools import partial

import numpy as np
import pandas as pd

from .. import utils
from .bam import read_genes
from .index import read_index, split_index

logger = logging.getLogger(__name__)

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
    (?P<T>[^,]*)\n
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


def read_counts(counts_path):
    """Read counts CSV as a pandas dataframe.

    :param counts_path: path to CSV
    :type counts_path: str

    :return: counts dataframe
    :rtype: pandas.DataFrame
    """
    dtypes = {'barcode': 'string', 'GX': 'string', **{column: np.uint8 for column in COLUMNS}}
    return pd.read_csv(counts_path, dtype=dtypes)


def read_counts_complemented(counts_path, genes_path):
    df_counts = read_counts(counts_path).merge(read_genes(genes_path)[['GX', 'strand']], on='GX')

    columns = ['barcode', 'GX'] + COLUMNS
    df_forward = df_counts[df_counts.strand == '+'][columns]
    df_reverse = df_counts[df_counts.strand == '-'][columns]

    df_reverse.columns = ['barcode', 'GX'] + CONVERSION_COLUMNS[::-1] + BASE_COLUMNS[::-1]
    df_reverse = df_reverse[columns]

    return pd.concat((df_forward, df_reverse))


def split_counts_by_umi(adata, df_counts, conversion='TC', filter_dict=None):
    filter_dict = filter_dict or {}
    for column, values in filter_dict.items():
        df_counts = df_counts[df_counts[column].isin(values)]

    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_id', inplace=True)

    new = np.zeros(adata.shape, dtype=int)
    counts = counts = df_counts[df_counts[conversion] > 0].groupby(['barcode', 'GX']).size()
    for (barcode, gene_id), count in counts.items():
        i = adata.obs.index.get_loc(barcode)
        j = adata.var.index.get_loc(gene_id)
        new[i, j] = count
    adata.layers['new_umi'] = new
    adata.layers['old_umi'] = adata.X - new

    adata.obs.reset_index(inplace=True)
    adata.var.reset_index(inplace=True)
    return adata


def count_conversions_part(
    conversions_path,
    counter,
    lock,
    pos,
    n_lines,
    snps=None,
    group_by=None,
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

        if group_by is None:
            return int(g['genome_i']) in snps.get(g['contig'], set())
        else:
            # TODO
            raise Exception()

    count_path = utils.mkstemp(dir=temp_dir)

    counts = None
    read_id = None
    count_base = True
    with open(conversions_path, 'r') as f, open(count_path, 'w') as out:
        f.seek(pos)

        groups = None
        prev_groups = None
        for i in range(n_lines):
            line = f.readline()
            prev_groups = groups
            groups = CONVERSIONS_PARSER.match(line).groupdict()

            if read_id != groups['read_id']:
                if read_id is not None and any(c > 0 for c in counts[:len(CONVERSION_IDX)]):
                    out.write(
                        f'{prev_groups["barcode"]},{prev_groups["umi"]},{prev_groups["GX"]},'
                        f'{",".join(str(c) for c in counts)}\n'
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

            if (i + 1) % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

        # Add last record
        if read_id is not None and any(c > 0 for c in counts[:len(CONVERSION_IDX)]):
            out.write(f'{groups["barcode"]},{groups["umi"]},{groups["GX"]},' f'{",".join(str(c) for c in counts)}\n')

    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return count_path


def count_conversions(
    conversions_path, index_path, counts_path, snps=None, group_by=None, quality=27, n_threads=8, temp_dir=None
):
    """Count the number of conversions of each read per barcode and gene, along with
    the total nucleotide content of the region each read mapped to, also per barcode.
    When a duplicate UMI for a barcode is observed, the read with the greatest
    number of conversions is selected.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param barcodes_path: path to write barcodes CSV
    :type barcodes_path: str
    :param genes_path: path to write genes CSV
    :type genes_path: str
    :param counts_path: path to write counts CSV
    :param counts_path: str
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional

    :return: (`barcodes_path`, `genes_path`, `counts_path`)
    :rtype: tuple
    """
    # Load index
    logger.debug(f'Loading index {index_path} for {conversions_path}')
    index = read_index(index_path)

    # Split index into n contiguous pieces
    logger.debug(f'Splitting index into {n_threads} parts')
    parts = split_index(index, n=n_threads)

    # Parse each part in a different process
    logger.debug(f'Spawning {n_threads} processes')
    n_lines = sum(idx[1] for idx in index)
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            count_conversions_part,
            conversions_path,
            counter,
            lock,
            snps=snps,
            group_by=group_by,
            quality=quality,
            temp_dir=tempfile.mkdtemp(dir=temp_dir)
        ), parts
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(async_result, counter, n_lines)
    pool.join()

    # Combine csvs
    combined_path = utils.mkstemp(dir=temp_dir)
    logger.debug(f'Combining intermediate parts to {combined_path}')
    with open(combined_path, 'wb') as out:
        for counts_part_path in async_result.get():
            with open(counts_part_path, 'rb') as f:
                shutil.copyfileobj(f, out)

    # Read in as dataframe and deduplicate based on CR and UB
    # If there are duplicates, select the row with the most conversions
    logger.debug(f'Deduplicating reads based on barcode and UMI to {counts_path}')
    df = pd.read_csv(
        combined_path,
        names=['barcode', 'umi', 'GX'] + COLUMNS,
        dtype={
            'barcode': 'string',
            'umi': 'string',
            'GX': 'string',
            **{column: np.uint8
               for column in COLUMNS}
        }
    )
    df_sorted = df.iloc[df[CONVERSION_COLUMNS].sum(axis=1).argsort()]
    df_filtered = df_sorted[~df_sorted.duplicated(subset=['barcode', 'umi'], keep='last')].drop(
        columns='umi'
    ).sort_values('barcode').reset_index(drop=True)
    df_filtered.to_csv(counts_path, index=False)

    return counts_path
