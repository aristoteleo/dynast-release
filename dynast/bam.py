import gzip
import logging
import os
import pickle
import time
import multiprocessing
import shutil
import tempfile
from collections import Counter
from functools import partial

import pysam
from tqdm import tqdm

from . import utils

logger = logging.getLogger(__name__)


def parse_read_contig(
    bam_path,
    counter,
    lock,
    contig,
    temp_dir=None,
    update_every=10000,
    process_every=10000,
):
    conversions_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)
    coverage_path = utils.mkstemp(dir=temp_dir)

    index = []

    keys_added = []
    coverage = {}
    mismatches = set()

    # NOTE: how to deal with duplicate UMIs?
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(conversions_path, 'w') as conversions_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = conversions_out.tell()
            # Add every 10000 to processed counter (updating every
            # iteration is very slow)
            # Can't use enumerate because n needs to be access outside from
            # the for loop.
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            # Skip reads that did not map to a gene
            if any(not read.has_tag(tag) for tag in ['GX', 'CB', 'UB']):
                continue

            read_id = read.query_name
            cr, ur = read.get_tag('CR'), read.get_tag('UR')
            cb, ub = read.get_tag('CB'), read.get_tag('UB')
            gx = read.get_tag('GX')
            strand = '-' if read.is_reverse else '+'
            sequence = read.seq.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()
            counts = Counter(reference)

            for read_i, genome_i, _genome_base in read.get_aligned_pairs(matches_only=False, with_seq=True):
                if genome_i is None:
                    continue

                key = (cr, genome_i)
                coverage.setdefault(key, 0)
                coverage[key] += 1
                keys_added.append(key)

                if read_i is None:
                    continue

                read_base = sequence[read_i]
                genome_base = _genome_base.upper()
                if 'N' in (genome_base, read_base):
                    continue

                if genome_base != read_base:
                    n_lines += 1
                    mismatches.add(key)
                    conversions_out.write(
                        f'{read_id},{cr},{cb},{ur},{ub},{gx},{strand},{contig},{genome_i},{genome_base},'
                        f'{read_base},{qualities[read_i]},{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )

            if n_lines > 0:
                index.append((pos, n_lines))

            # Remove coverage without any mismatches.
            # We do this to save memory.
            # We also don't do this every iteration in interest of runtime.
            if n % process_every == 0:
                new_keys_added = []
                for key in list(set(keys_added)):
                    if key[1] < read.reference_start:
                        if key not in mismatches:
                            del coverage[key]
                    else:
                        new_keys_added.append(key)
                mismatches = mismatches.intersection(new_keys_added)
                keys_added = new_keys_added

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    # Save index
    with gzip.open(index_path, 'wb') as f:
        pickle.dump(index, f)

    # Clean coverage before writing
    for key in list(set(keys_added)):
        if key not in mismatches:
            del coverage[key]
    with open(coverage_path, 'w') as f:
        for (cr, genome_i), count in coverage.items():
            f.write(f'{cr},{contig},{genome_i},{count}\n')

    return conversions_path, index_path, coverage_path


def parse_all_reads(bam_path, conversions_path, index_path, coverage_path, n_threads=8, temp_dir=None):
    """
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        contigs = bam.references
        n_reads = sum(index.total for index in bam.get_index_statistics())

    # Initialize and run pool
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.map_async(
        partial(parse_read_contig, bam_path, counter, lock, temp_dir=tempfile.mkdtemp(dir=temp_dir)), contigs
    )

    # Display progres bar
    pbar = tqdm(total=n_reads, ascii=True)
    previous_progress = 0
    while not async_result.ready():
        time.sleep(0.05)
        progress = counter.value
        pbar.update(progress - previous_progress)
        previous_progress = progress
    pbar.close()

    # Combine csvs
    pos = 0
    index = []
    with open(conversions_path, 'wb') as conversions_out, \
        open(coverage_path, 'wb') as coverage_out:
        for conversions_part_path, index_part_path, coverage_part_path in async_result.get():
            with open(conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, conversions_out)
            with open(coverage_part_path, 'rb') as f:
                shutil.copyfileobj(f, coverage_out)

            # Parse index
            with gzip.open(index_part_path, 'rb') as f:
                index_part = pickle.load(f)
            for p, n in index_part:
                index.append((pos + p, n))
            pos += os.path.getsize(conversions_part_path)
    with gzip.open(index_path, 'wb') as index_out:
        pickle.dump(index, index_out)

    return conversions_path, index_path, coverage_path
