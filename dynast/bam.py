import os
import time
import multiprocessing
import shutil
import tempfile
from collections import Counter
from functools import partial

import pysam
from tqdm import tqdm

BASE_IDX = {
    'A': 0,
    'C': 1,
    'G': 2,
    'T': 3,
}
IDX_BASE = {v: k for k, v in BASE_IDX.items()}
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
IDX_CONVERSION = {v: k for k, v in CONVERSION_IDX.items()}


def parse_read_contig(
    bam_path, processed, lock, contig, temp_dir=None, update_every=10000, process_every=100000,
):
    conversions_dir = os.path.join(temp_dir, 'conversions')
    info_dir = os.path.join(temp_dir, 'info')
    coverage_dir = os.path.join(temp_dir, 'coverage')
    os.makedirs(conversions_dir, exist_ok=True)
    os.makedirs(info_dir, exist_ok=True)
    os.makedirs(coverage_dir, exist_ok=True)

    conversions_path = os.path.join(conversions_dir, f'{contig}.csv')
    info_path = os.path.join(info_dir, f'{contig}.csv')
    coverage_path = os.path.join(coverage_dir, f'{contig}.csv')

    keys_added = []
    coverage = {}
    mismatches = set()

    # NOTE: how to deal with duplicate UMIs?
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(info_path, 'w') as info_out, \
        open(conversions_path, 'w') as conversions_out:
        for read in bam.fetch(contig):
            # Add every 10000 to processed counter (updating every
            # iteration is very slow)
            # Can't use enumerate because n needs to be access outside from
            # the for loop.
            n += 1
            if n % update_every == 0:
                lock.acquire()
                processed.value += update_every
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

            mismatched = False
            for read_i, genome_i, _genome_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                read_base = sequence[read_i]
                genome_base = _genome_base.upper()
                if 'N' in (genome_base, read_base):
                    continue

                key = (cr, genome_i)
                coverage.setdefault(key, 0)
                coverage[key] += 1
                keys_added.append(key)

                if genome_base != read_base:
                    mismatched = True
                    mismatches.add(key)
                    conversions_out.write(
                        f'{read_id},{contig},{strand},{genome_i},{genome_base},{read_base},{qualities[read_i]}\n'
                    )

            if mismatched:
                counts = Counter(reference)
                info_out.write(f'{read_id},{cr},{cb},{ur},{ub},{gx},{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n')

            # Remove coverage without any mismatches.
            # We do this to save memory.
            # We also don't do this every iteration in interest of runtime.
            if n % process_every == 0:
                new_keys_added = []
                for key in keys_added:
                    if key[1] < read.reference_start:
                        if key not in mismatches and key in coverage:
                            del coverage[key]
                    else:
                        new_keys_added.append(key)
                mismatches = mismatches.intersection(new_keys_added)
                keys_added = new_keys_added

        lock.acquire()
        processed.value += n % update_every
        lock.release()

    # Clean coverage before writing
    for key in keys_added:
        if key not in mismatches and key in coverage:
            del coverage[key]
    with open(coverage_path, 'w') as f:
        for (cr, genome_i), count in coverage.items():
            f.write(f'{cr},{contig},{genome_i},{count}\n')

    return info_path, conversions_path, coverage_path

def parse_all_reads(bam_path, info_path, conversions_path, coverage_path, n_threads=8, temp_dir=None):
    """
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        contigs = bam.references

    # Initialize and run pool
    manager = multiprocessing.Manager()
    processed = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.map_async(
        partial(parse_read_contig, bam_path, processed, lock, temp_dir=tempfile.mkdtemp(dir=temp_dir)),
        contigs
    )

    # Display progres bar
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        pbar = tqdm(total=sum(index.total for index in bam.get_index_statistics()), ascii=True)
    previous_progress = 0
    while not async_result.ready():
        time.sleep(0.05)
        progress = processed.value
        pbar.update(progress - previous_progress)
        previous_progress = progress
    pbar.close()

    # Combine csvs
    with open(info_path, 'wb') as info_out, \
        open(conversions_path, 'wb') as conversions_out, \
        open(coverage_path, 'wb') as coverage_out:
        for info_part_path, conversions_part_path, coverage_part_path in async_result.get():
            with open(info_part_path, 'rb') as f:
                shutil.copyfileobj(f, info_out)
            with open(conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, conversions_out)
            with open(coverage_part_path, 'rb') as f:
                shutil.copyfileobj(f, coverage_out)

    return info_path, conversions_path, coverage_path
