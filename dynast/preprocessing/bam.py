import logging
import os
import shutil
import tempfile
from collections import Counter
from functools import partial

import numpy as np
import pandas as pd
import pysam

from .. import utils

logger = logging.getLogger(__name__)


def read_conversions(conversions_path):
    df = pd.read_csv(
        conversions_path,
        dtype={
            'contig': 'category',
            'genome_i': int
        },
        usecols=['contig', 'genome_i'],
    ).drop_duplicates()

    conversions = {}
    for contig, idx in df.groupby('contig').indices.items():
        conversions[contig] = set(df.iloc[idx]['genome_i'])
    return conversions


def read_no_conversions(no_conversions_path):
    df = pd.read_csv(
        no_conversions_path,
        dtype={
            'read_id': 'string',
            'barcode': 'string',
            'umi': 'string',
            'A': np.uint8,
            'C': np.uint8,
            'G': np.uint8,
            'T': np.uint8,
        }
    )
    return df


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


def parse_read_contig(
    bam_path,
    counter,
    lock,
    contig,
    barcodes=None,
    read_group_as_barcode=False,
    extract_genes=False,
    use_corrected=False,
    temp_dir=None,
    update_every=50000,
    nasc=False
):
    """Parse all reads mapped to a contig, outputing conversion
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

    :return: (`conversions_path`, `index_path`)
             `conversions_path`: path to temporary conversions CSV
             `index_path`: path to temporary index pickle
    :rtype: tuple
    """
    conversions_path = utils.mkstemp(dir=temp_dir)
    no_conversions_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)

    # Index will be used later for fast splitting
    # Will contain a tuple of (file position, number of lines) for every read
    # that has at least one mismatch.
    index = []
    genes = {}
    paired = {}

    required_tags = ['GX']
    if not read_group_as_barcode:
        required_tags.append('UB')
        if use_corrected:
            required_tags.append('CB')

    # Can't use enumerate because n needs to be access outside from
    # the for loop.
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam, \
        open(conversions_path, 'w') as conversions_out, \
        open(no_conversions_path, 'w') as no_conversions_out:
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
            if any(not read.has_tag(tag) for tag in required_tags) or read.is_duplicate:
                continue

            # Extract read and reference information
            read_id = read.query_name
            barcode = read.get_tag('RG') if read_group_as_barcode else (
                read.get_tag('CB') if use_corrected else read.get_tag('CR')
            )

            # If barcodes are provided, skip any barcodes not in the list.
            if barcodes and barcode not in barcodes:
                continue

            umi = None if read_group_as_barcode else read.get_tag('UB')
            gx = read.get_tag('GX') if read.has_tag('GX') else ''
            sequence = read.seq.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()

            if extract_genes:
                gn = read.get_tag('GN') if read.has_tag('GN') else ''
                strand = '-' if read.is_reverse else '+'
                genes.setdefault(gx, (gn, strand))

            # Count number of nucleotides in the region this read mapped to
            counts = Counter(reference)

            # Used exclusively for paired reads. So that we don't count conversions
            # in overlapping regions twice.
            # key: genome_i, value: (quality, original, converted)
            conversions = {}

            # If we are dealing with paired-end reads, we have to check for
            # overlaps.
            if read.is_paired:
                if read_id not in paired:
                    paired[read_id] = (
                        set(read.get_reference_positions()),
                        reference,
                        sequence,
                        qualities,
                        barcode,
                        umi,
                        gx,
                        read.get_aligned_pairs(matches_only=True, with_seq=True),
                    )
                    continue

                (
                    mate_reference_positions,
                    mate_reference,
                    mate_sequence,
                    mate_qualities,
                    mate_barcode,
                    mate_umi,
                    mate_gx,
                    aligned_pairs,
                ) = paired[read_id]

                # Update counts. union(A, B) = A + B - intersection(A, B)
                counts += Counter(mate_reference)
                read_reference_positions = read.get_reference_positions()
                overlap_pos = list(mate_reference_positions.intersection(read_reference_positions))
                overlap = Counter()
                for base, i in zip(reference, read_reference_positions):
                    if i in overlap_pos:
                        overlap[base] += 1
                counts -= overlap

                for read_i, genome_i, _genome_base in aligned_pairs:
                    read_base = mate_sequence[read_i]
                    genome_base = _genome_base.upper()
                    if 'N' in (genome_base, read_base):
                        continue

                    if genome_base != read_base and (not nasc or genome_i not in overlap_pos):
                        # Add this conversion to the conversions dictionary.
                        # We don't write this conversion immediately because
                        # we want to take the conversion with higher quality
                        # if conversions exist in same positions in overlapping
                        # segments.
                        conversions[genome_i] = (mate_qualities[read_i], genome_base, read_base)
                del paired[read_id]

            # Iterate through every mapped position.
            for read_i, genome_i, _genome_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
                read_base = sequence[read_i]
                genome_base = _genome_base.upper()
                if 'N' in (genome_base, read_base):
                    continue

                quality = qualities[read_i]
                # If this index exists in the conversions dictionary,
                # then select the conversion with higher quality
                if genome_i in conversions:
                    if not nasc:
                        mate_quality, _, mate_read_base = conversions[genome_i]
                        if mate_quality > quality:
                            read_base = mate_read_base
                            quality = mate_quality
                    del conversions[genome_i]

                if genome_base != read_base:
                    n_lines += 1
                    if read_group_as_barcode:
                        conversions_out.write(
                            f'{read_id},{barcode},{gx},{contig},{genome_i},'
                            f'{genome_base},{read_base},{quality},'
                            f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                        )
                    else:
                        conversions_out.write(
                            f'{read_id},{barcode},{umi},{gx},{contig},{genome_i},'
                            f'{genome_base},{read_base},{quality},'
                            f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                        )

            # Flush cached conversions
            for genome_i, (quality, genome_base, read_base) in conversions.items():
                n_lines += 1
                if read_group_as_barcode:
                    conversions_out.write(
                        f'{read_id},{barcode},{gx},{contig},{genome_i},'
                        f'{genome_base},{read_base},{quality},'
                        f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )
                else:
                    conversions_out.write(
                        f'{read_id},{barcode},{umi},{gx},{contig},{genome_i},'
                        f'{genome_base},{read_base},{quality},'
                        f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )

            # Add to index if lines were written
            if n_lines > 0:
                index.append((pos, n_lines))
            else:
                # Otherwise, this read did not contain any conversions.
                if read_group_as_barcode:
                    no_conversions_out.write(
                        f'{read_id},{barcode},{gx},{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )
                else:
                    no_conversions_out.write(
                        f'{read_id},{barcode},{umi},{gx},{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]}\n'
                    )

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    # Save index
    index_path = utils.write_pickle(index, index_path, protocol=4)

    return (conversions_path, index_path, no_conversions_path,
            genes) if extract_genes else (conversions_path, index_path, no_conversions_path)


def parse_all_reads(
    bam_path,
    conversions_path,
    index_path,
    no_conversions_path,
    barcodes=None,
    genes_path=None,
    read_group_as_barcode=False,
    use_corrected=False,
    n_threads=8,
    temp_dir=None,
    nasc=False,
):
    """Parse all reads in the provided BAM file. Read conversion and coverage
    information is outputed to the provided corresponding paths.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions_path: path to CSV to write conversion information
    :type conversions_path: str
    :param index_path: path to index pickle to write
    :type index_path: str
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional

    :return: (`conversions_path`, `index_path`)
    :rtype: tuple
    """
    logger.debug(f'Extracting contigs from BAM {bam_path}')
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    # Initialize and run pool
    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.map_async(
        partial(
            parse_read_contig,
            bam_path,
            counter,
            lock,
            barcodes=barcodes,
            read_group_as_barcode=read_group_as_barcode,
            extract_genes=genes_path is not None,
            use_corrected=use_corrected,
            temp_dir=tempfile.mkdtemp(dir=temp_dir),
            nasc=nasc
        ), contigs
    )
    pool.close()

    # Display progres bar
    logger.debug(f'Processing contigs {contigs} from BAM')
    utils.display_progress_with_counter(async_result, counter, n_reads)
    pool.join()

    # Combine csvs
    logger.debug(f'Writing conversions to {conversions_path}')
    pos = 0
    index = []
    genes = {}
    with open(conversions_path, 'wb') as conversions_out, \
        open(no_conversions_path, 'wb') as no_conversions_out:
        if read_group_as_barcode:
            conversions_out.write(b'read_id,barcode,GX,contig,genome_i,original,converted,quality,A,C,G,T\n')
            no_conversions_out.write(b'read_id,barcode,GX,A,C,G,T\n')
        else:
            conversions_out.write(b'read_id,barcode,umi,GX,contig,genome_i,original,converted,quality,A,C,G,T\n')
            no_conversions_out.write(b'read_id,barcode,umi,GX,A,C,G,T\n')

        pos = conversions_out.tell()

        for parts in async_result.get():
            conversions_part_path = parts[0]
            index_part_path = parts[1]
            no_conversions_part_path = parts[2]
            if genes_path is not None:
                genes_part = parts[3]
                genes.update(genes_part)
            with open(conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, conversions_out)
            with open(no_conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, no_conversions_out)

            # Parse index
            index_part = utils.read_pickle(index_part_path)
            for p, n in index_part:
                index.append((pos + p, n))
            pos += os.path.getsize(conversions_part_path)
    logger.debug(f'Writing conversions index to {index_path}')
    index_path = utils.write_pickle(index, index_path, protocol=4)

    if genes_path is not None:
        logger.debug(f'Writing genes CSV to {genes_path}')
        with open(genes_path, 'w') as f:
            f.write('GX,GN,strand\n')
            for gx in sorted(genes.keys()):
                gn, strand = genes[gx]
                f.write(f'{gx},{gn},{strand}\n')

    return (conversions_path, index_path, no_conversions_path,
            genes_path) if genes_path is not None else (conversions_path, index_path, no_conversions_path)
