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

from .. import utils

logger = logging.getLogger(__name__)


def parse_read_contig(
    bam_path,
    counter,
    lock,
    contig,
    use_corrected=False,
    temp_dir=None,
    update_every=5000,
    process_every=10000,
):
    """Parse all reads mapped to a contig, outputing conversion and coverage
    information as temporary CSVs. This function is designed to be called as a
    separate process.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param contig: only reads that map to this contig will be processed
    :type contig: str
    :param use_corrected: whether or not to require the corrected barcode (CB)
                          tag to be present. If `False`, reads without this tag
                          will be skipped, defaults to `False`
    :type use_corrected: bool, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `5000`
    :type update_every: int, optional
    :param process_every: update coverage information every this many reads.
                          This is in place so that unneeded coverage information
                          is flushed from the buffer, in the interest of memory.
                          Raising this value will result in faster runtimes but
                          higher memory usage. A value of `10000` requires 16 GB
                          of memory on average, defaults to `10000`
    :type process_every: int, optional

    :return: (`conversions_path`, `index_path`, `coverage_path`)
             `conversions_path`: path to temporary conversions CSV
             `index_path`: path to temporary index pickle
             `converage_path`: path to temporary coverage CSV
    :rtype: tuple
    """
    conversions_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)
    coverage_path = utils.mkstemp(dir=temp_dir)

    # Index will be used later for fast splitting
    # Will contain a tuple of (file position, number of lines) for every read
    # that has at least one mismatch.
    index = []

    # Keep track of keys and mismatches added since last time they were processed
    keys_added = []
    mismatches = set()

    coverage = {}
    required_tags = ['GX', 'UB']
    if use_corrected:
        required_tags.append('CB')

    # Can't use enumerate because n needs to be access outside from
    # the for loop.
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(conversions_path, 'w') as conversions_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = conversions_out.tell()
            # Update every some interval. Updating every loop is very slow.
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            # Skip reads that do not contain these tags
            if any(not read.has_tag(tag) for tag in required_tags):
                continue

            # Extract read and reference information
            read_id = read.query_name
            cr, ur = read.get_tag('CR'), read.get_tag('UR')
            cb = read.get_tag('CB')
            ub = read.get_tag('UB') if read.has_tag('UB') else ''
            gx = read.get_tag('GX')
            gn = read.get_tag('GN') if read.has_tag('GN') else ''
            strand = '-' if read.is_reverse else '+'
            sequence = read.seq.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()

            # Count number of nucleotides in the region this read mapped to
            counts = Counter(reference)

            # Iterate through every mapped position.
            # Note that matches_only=False because we need to extract coverage
            # information, and in order to do that we need to consider gaps as well.
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
                        f'{read_id},{cr},{cb},{ur},{ub},{gx},{gn},{strand},{contig},{genome_i},{genome_base},'
                        f'{read_base},{qualities[read_i]},{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )

            # Add to index if lines were written
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
        pickle.dump(index, f, protocol=4)

    # Clean coverage before writing
    for key in list(set(keys_added)):
        if key not in mismatches:
            del coverage[key]
    with open(coverage_path, 'w') as f:
        for (cr, genome_i), count in coverage.items():
            f.write(f'{cr},{contig},{genome_i},{count}\n')

    return conversions_path, index_path, coverage_path


def parse_all_reads(bam_path, conversions_path, index_path, coverage_path, n_threads=8, temp_dir=None):
    """Parse all reads in the provided BAM file. Read conversion and coverage
    information is outputed to the provided corresponding paths.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions_path: path to CSV to write conversion information
    :type conversions_path: str
    :param index_path: path to index pickle to write
    :type index_path: str
    :param coverage_path: path to CSV to write coverage information
    :type coverage_path: str
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional

    :return: (`conversions_path`, `index_path`, `coverage_path`)
    :rtype: tuple
    """
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    # Initialize and run pool
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.map_async(
        partial(parse_read_contig, bam_path, counter, lock, temp_dir=tempfile.mkdtemp(dir=temp_dir)), contigs
    )
    pool.close()

    # Display progres bar
    with tqdm(unit='reads', total=n_reads, ascii=True, unit_scale=True) as pbar:
        previous_progress = 0
        while not async_result.ready():
            time.sleep(0.05)
            progress = counter.value
            pbar.update(progress - previous_progress)
            previous_progress = progress
        pool.join()

    # Combine csvs
    pos = 0
    index = []
    with open(conversions_path, 'wb') as conversions_out, \
        open(coverage_path, 'wb') as coverage_out:
        conversions_out.write(b'read_id,CR,CB,UR,UB,GX,GN,strand,contig,genome_i,original,converted,A,C,G,T\n')
        coverage_out.write(b'CR,contig,genome_i,coverage\n')
        pos = conversions_out.tell()

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
        pickle.dump(index, index_out, protocol=4)

    return conversions_path, index_path, coverage_path
