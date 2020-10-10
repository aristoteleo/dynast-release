import gzip
import logging
import os
import multiprocessing
import pickle
import re
import shutil
import tempfile
import time
from functools import partial

from tqdm import tqdm

from . import utils

logger = logging.getLogger(__name__)

PARSER = re.compile(r'''^
    (?P<read_id>[^,]*),
    (?P<CR>[^,]*),
    (?P<CB>[^,]*),
    (?P<UR>[^,]*),
    (?P<UB>[^,]*),
    (?P<GX>[^,]*),
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
$''', re.VERBOSE)

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

def split_index(index, n=8):
    """
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

def count_conversions_part(conversions_path, counter, lock, pos, n_lines, quality=27, temp_dir=None, update_every=10000):
    """
    """
    count_path = utils.mkstemp(dir=temp_dir)

    counts = None
    read_id = None
    count_base = True
    with open(conversions_path, 'r') as f, open(count_path, 'w') as out:
        f.seek(pos)

        for i in range(n_lines):
            line = f.readline()
            groups = PARSER.match(line).groupdict()

            if read_id != groups['read_id']:
                if read_id is not None and any(c > 0 for c in counts[:len(CONVERSION_IDX)]):
                    out.write(
                        f'{read_id},{groups["CR"]},{groups["GX"]},{groups["strand"]},{groups["UB"]},'
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
    lock.acquire()
    counter.value += n_lines % update_every
    lock.release()

    return count_path


def count_conversions(conversions_path, index_path, count_path, quality=27, n_threads=8, temp_dir=None):
    """
    """
    # Load index
    with gzip.open(index_path, 'rb') as f:
        index = pickle.load(f)

    # Split index into n contiguous pieces
    parts = split_index(index, n=n_threads)

    # Parse each part in a different process
    n_lines = sum(idx[1] for idx in index)
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.starmap_async(
        partial(count_conversions_part, conversions_path, counter, lock, quality=quality, temp_dir=tempfile.mkdtemp(dir=temp_dir)),
        parts
    )

    # Display progres bar
    pbar = tqdm(total=n_lines, ascii=True)
    previous_progress = 0
    while not async_result.ready():
        time.sleep(0.05)
        progress = counter.value
        pbar.update(progress - previous_progress)
        previous_progress = progress
    pbar.close()

    # Combine csvs
    with open(count_path, 'wb') as out:
        for count_part_path in async_result.get():
            with open(count_part_path, 'rb') as f:
                shutil.copyfileobj(f, out)

    return count_path
