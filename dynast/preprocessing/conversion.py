import gzip
import logging
import multiprocessing
import pickle
import re
import shutil
import tempfile
import time
from functools import partial

import numpy as np
import pandas as pd
from tqdm import tqdm

from .. import utils

logger = logging.getLogger(__name__)

CONVERSIONS_PARSER = re.compile(
    r'''^
    (?P<read_id>[^,]*),
    (?P<CR>[^,]*),
    (?P<CB>[^,]*),
    (?P<UR>[^,]*),
    (?P<UB>[^,]*),
    (?P<GX>[^,]*),
    (?P<GN>[^,]*),
    (?P<strand>[^,]*),
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


def read_genes(genes_path):
    """Read genes CSV as a pandas dataframe.

    :param genes_path: path to CSV
    :type genes_path: str

    :return: genes dataframe
    :rtype: pandas.DataFrame
    """
    dtypes = {
        'GX': 'string',
        'GN': 'string',
        'strand': 'category',
    }
    return pd.read_csv(genes_path, dtype=dtypes)


def split_index(index, n=8):
    """Split a conversions index, which is a list of tuples (file position,
    number of lines), one for each read, into `n` approximately equal parts.
    This function is used to split the conversions CSV for multiprocessing.

    :param index: conversions index
    :type index: list
    :param n: number of splits, defaults to `8`
    :type n: int, optional

    :return: list of parts, where each part is a (file position, number of lines) tuple
    :rtype: list
    """
    n_lines = sum(idx[1] for idx in index)
    target = (n_lines // n) + 1  # add one to prevent underflow

    # Split the index to "approximately" equal parts
    parts = []
    start_pos = None
    current_size = 0
    for pos, size in index:
        if start_pos is None:
            start_pos = pos
        current_size += size

        if current_size >= target:
            parts.append((start_pos, current_size))
            start_pos = None
            current_size = 0
    if current_size > 0:
        parts.append((start_pos, current_size))

    return parts


def count_conversions_part(
    conversions_path, counter, lock, pos, n_lines, use_corrected=False, quality=27, temp_dir=None, update_every=10000
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
    :param use_corrected: whether or not to use corrected barcodes in the `barcode`
                          column of the counts CSV, defaults to `False`
    :type use_corrected: bool, optional
    :param quality: only count conversions with PHRED quality greater than this value,
                    defaults to `27`
    :type quality: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `10000`
    :type update_every: int, optional

    :return: (`count_path`, `barcodes`, `genes`)
             `count_path`: path to temporary counts CSV
             `barcodes`: dictionary of raw to corrected barcode mappings
             `genes`: dictionary of gene ID to (gene name, strand) mappings
    :rtype: tuple
    """
    count_path = utils.mkstemp(dir=temp_dir)

    barcodes = {}
    genes = {}
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
                    barcodes.setdefault(prev_groups["CR"], prev_groups["CB"])
                    genes.setdefault(prev_groups["GX"], (prev_groups['GN'], prev_groups['strand']))

                    barcode = prev_groups["CB"] if use_corrected else prev_groups["CR"]
                    out.write(
                        f'{barcode},{prev_groups["UB"]},{prev_groups["GX"]},'
                        f'{",".join(str(c) for c in counts)}\n'
                    )
                counts = [0] * (len(CONVERSION_IDX) + len(BASE_IDX))
                read_id = groups['read_id']
                count_base = True
            if int(groups['quality']) > quality:
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
            barcodes.setdefault(groups["CR"], groups["CB"])
            genes.setdefault(groups["GX"], (groups['GN'], groups['strand']))

            barcode = groups["CB"] if use_corrected else groups["CR"]
            out.write(f'{barcode},{groups["UB"]},{groups["GX"]},' f'{",".join(str(c) for c in counts)}\n')

    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return count_path, barcodes, genes


def count_conversions(
    conversions_path,
    index_path,
    barcodes_path,
    genes_path,
    counts_path,
    use_corrected=False,
    quality=27,
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
    :param barcodes_path: path to write barcodes CSV
    :type barcodes_path: str
    :param genes_path: path to write genes CSV
    :type genes_path: str
    :param counts_path: path to write counts CSV
    :param counts_path: str
    :param use_corrected: whether or not to use corrected barcodes in the `barcode`
                          column of the counts CSV, defaults to `False`
    :type use_corrected: bool, optional
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
    with gzip.open(index_path, 'rb') as f:
        index = pickle.load(f)

    # Split index into n contiguous pieces
    logger.debug(f'Splitting index into {n_threads} parts')
    parts = split_index(index, n=n_threads)

    # Parse each part in a different process
    logger.debug(f'Spawning {n_threads} processes')
    n_lines = sum(idx[1] for idx in index)
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.starmap_async(
        partial(
            count_conversions_part,
            conversions_path,
            counter,
            lock,
            use_corrected=use_corrected,
            quality=quality,
            temp_dir=tempfile.mkdtemp(dir=temp_dir)
        ), parts
    )
    pool.close()

    # Display progres bar
    with tqdm(unit='reads', total=n_lines, ascii=True, unit_scale=True) as pbar:
        previous_progress = 0
        while not async_result.ready():
            time.sleep(0.05)
            progress = counter.value
            pbar.update(progress - previous_progress)
            previous_progress = progress
    pool.join()

    # Combine csvs
    barcodes = {}
    genes = {}
    combined_path = utils.mkstemp(dir=temp_dir)
    logger.debug(f'Combining intermediate parts to {combined_path}')
    with open(combined_path, 'wb') as out:
        for counts_part_path, barcodes_part, genes_part in async_result.get():
            barcodes.update(barcodes_part)
            genes.update(genes_part)
            with open(counts_part_path, 'rb') as f:
                shutil.copyfileobj(f, out)

    # Write barcode and gene mapping
    logger.debug(f'Writing barcodes to {barcodes_path}')
    with open(barcodes_path, 'w') as f:
        f.write('CR,CB\n')
        for cr in sorted(barcodes.keys()):
            f.write(f'{cr},{barcodes[cr]}\n')
    logger.debug(f'Writing genes to {genes_path}')
    with open(genes_path, 'w') as f:
        f.write('GX,GN,strand\n')
        for gx in sorted(genes.keys()):
            gn, strand = genes[gx]
            f.write(f'{gx},{gn},{strand}\n')

    # Read in as dataframe and deduplicate based on CR and UB
    # If there are duplicates, select the row with the most conversions
    logger.debug(f'Deduplicating reads based on barcode and UMI to {counts_path}')
    df = pd.read_csv(
        combined_path,
        names=['barcode', 'UB', 'GX'] + COLUMNS,
        dtype={
            'barcode': 'string',
            'CB': 'string',
            'GX': 'string',
            **{column: np.uint8
               for column in COLUMNS}
        }
    )
    df_sorted = df.iloc[df[CONVERSION_COLUMNS].sum(axis=1).argsort()]
    df_filtered = df_sorted[~df_sorted.duplicated(subset=['barcode', 'UB'], keep='last')].drop(
        columns='UB'
    ).sort_values('barcode').reset_index(drop=True)
    df_filtered.to_csv(counts_path, index=False)

    return barcodes_path, genes_path, counts_path
