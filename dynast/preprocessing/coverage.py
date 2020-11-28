import gzip
import logging
import os
import pickle
import time
import multiprocessing
import shutil
import tempfile
from functools import partial

import pysam
from tqdm import tqdm

from .. import utils

logger = logging.getLogger(__name__)


def calculate_coverage_contig(
    bam_path,
    counter,
    lock,
    contig,
    indices,
    barcodes=None,
    read_group_as_barcode=False,
    use_corrected=False,
    temp_dir=None,
    update_every=50000,
    flush_every=10000
):
    coverage_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)

    # Index will be used for fast splitting, just like for conversions.csv
    index = []

    coverage = {}

    required_tags = ['GX']
    if not read_group_as_barcode:
        required_tags.append('UB')
    if use_corrected:
        required_tags.append('CB')

    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(coverage_path, 'w') as coverage_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = coverage_out.tell()
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

            if any(not read.has_tag(tag) for tag in required_tags):
                continue

            barcode = read.get_tag('RG') if read_group_as_barcode else (
                read.get_tag('CB') if use_corrected else read.get_tag('CR')
            )
            if barcodes and barcode not in barcodes:
                continue

            for genome_i in range(read.reference_start, read.reference_end):
                if genome_i in indices:
                    key = (barcode, genome_i)
                    coverage.setdefault(key, 0)
                    coverage[key] += 1

            if n % flush_every == 0:
                for key in sorted(list(coverage.keys())):
                    if key[1] < read.reference_start:
                        n_lines += 1
                        coverage_out.write(f'{key[0]},{contig},{key[1]},{coverage[key]}\n')
                        del coverage[key]

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
    barcodes=None,
    read_group_as_barcode=False,
    use_corrected=False,
    n_threads=8,
    temp_dir=None,
):
    logger.debug(f'Extracting contigs from BAM {bam_path}')
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    logger.debug(f'Spawning {n_threads} processes')
    manager = multiprocessing.Manager()
    counter = manager.Value('I', 0)
    lock = manager.Lock()
    pool = multiprocessing.Pool(n_threads)
    async_result = pool.starmap_async(
        partial(
            calculate_coverage_contig,
            bam_path,
            counter,
            lock,
            barcodes=barcodes,
            read_group_as_barcode=read_group_as_barcode,
            use_corrected=use_corrected,
            temp_dir=tempfile.mkdtemp(dir=temp_dir)
        ), [(contig, conversions.get(contig, set())) for contig in contigs]
    )
    pool.close()

    # Display progres bar
    logger.debug(f'Processing contigs {contigs} from BAM')
    with tqdm(unit='reads', total=n_reads, ascii=True, unit_scale=True) as pbar:
        previous_progress = 0
        while not async_result.ready():
            time.sleep(0.05)
            progress = counter.value
            pbar.update(progress - previous_progress)
            previous_progress = progress
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
