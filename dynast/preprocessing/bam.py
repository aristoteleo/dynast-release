import logging
import os
import shutil
import tempfile
from collections import Counter
from functools import partial

import numpy as np
import pandas as pd
import pysam

from . import gtf
from .. import config, utils

logger = logging.getLogger(__name__)


def read_conversions(conversions_path, *args, **kwargs):
    """Read conversions CSV as a pandas DataFrame.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    :param conversions_path: path to conversions CSV
    :type conversions_path: str

    :return: conversions dataframe
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(
        conversions_path,
        dtype={
            'read_id': 'string',
            'barcode': 'string',
            'umi': 'string',
            'contig': 'category',
            'genome_i': np.uint32,
            'original': 'string',
            'converted': 'string',
            'A': np.uint8,
            'C': np.uint8,
            'G': np.uint8,
            'T': np.uint8,
            'velocity': 'category',
            'transcriptome': bool,
        },
        *args,
        **kwargs
    )
    return df


def parse_read_contig(
    bam_path,
    counter,
    lock,
    contig,
    gene_infos=None,
    transcript_infos=None,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    barcodes=None,
    temp_dir=None,
    update_every=2000,
    nasc=False,
    velocity=True
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
    :param gene_infos: dictionary containing gene information, as returned by
                       `preprocessing.gtf.parse_gtf`, required if `velocity=True`,
                       defaults to `None`
    :type gene_infos: dictionary
    :param transcript_infos: dictionary containing transcript information,
                             as returned by `preprocessing.gtf.parse_gtf`,
                             required if `velocity=True`, defaults to `None`
    :type transcript_infos: dictionary
    :param strand: strandedness of the sequencing protocol, defaults to `forward`,
                   may be one of the following: `forward`, `reverse`, `None` (unstranded)
    :type strand: str, optional
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param update_every: update the counter every this many reads, defaults to `5000`
    :type update_every: int, optional
    :param nasc: flag to change behavior to match NASC-seq pipeline, defaults to `False`
    :type nasc: bool, optional
    :param velocity: whether or not to assign a velocity type to each read,
                     defaults to `True`
    :type velocity: bool, optional

    :return: (path to conversions, path to conversions index,
              path to no conversions, path to no conversions index)
    :rtype: (str, str, str, str)
    """
    # Sort genes by segment positions
    if velocity:
        gene_order = sorted(gene_infos.keys(), key=lambda gene: tuple(gene_infos[gene]['segment']))

    def assign_velocity_type(positions, gx=None, read_strand=None):
        alignment = gtf.SegmentCollection.from_positions(positions)

        # STAR assigns a GX tag to a read iff it only overlaps exons of one transcript model
        # But this doesn't mean it's still true for other transcript models of the same gene
        # If the alignment overlaps introns of any transcript model of the given gene,
        # it is ambiguous
        if gx:
            for transcript_id in gene_infos[gx]['transcripts']:
                if alignment.is_overlapping(transcript_infos[transcript_id]['introns']):
                    return gx, 'ambiguous'
            return gx, 'spliced'
        else:
            # Read was not assigned a GX tag, which means it overlaps intronic
            # regions, but we need to find which gene it is
            # note that the read could also map to an intergenic region, in which
            # case we ignore the read
            assigned_gene = None
            assigned_type = None
            for gene in gene_order:
                # On different strand
                if read_strand and read_strand != gene_infos[gene]['strand']:
                    continue
                # Out of range
                gene_segment = gene_infos[gene]['segment']
                if alignment.start >= gene_segment.end:
                    continue
                elif alignment.end <= gene_segment.start:
                    break

                any_exon_only = False
                any_intron_overlap = False
                for transcript_id in gene_infos[gene]['transcripts']:
                    transcript = transcript_infos[transcript_id]
                    # Check only if we have to
                    if not any_exon_only:
                        any_exon_only |= alignment.is_subset(transcript['exons'])
                    if not any_intron_overlap:
                        any_intron_overlap |= alignment.is_overlapping(transcript['introns'])

                if assigned_gene and (any_exon_only or any_intron_overlap):
                    return None, None

                # These consider every possible combination of
                # any_exon_only and any_intron_overlap
                if not any_exon_only and not any_intron_overlap:
                    continue
                if any_exon_only and not any_intron_overlap:
                    assigned_gene, assigned_type = gene, 'spliced'
                elif not any_exon_only and any_intron_overlap:
                    assigned_gene, assigned_type = gene, 'unspliced'
                else:  # any_exon_only = any_intron_overlap = True
                    assigned_gene, assigned_type = gene, 'ambiguous'

            return assigned_gene, assigned_type

        # SHOULD NEVER REACH HERE
        raise Exception()

    conversions_path = utils.mkstemp(dir=temp_dir)
    no_conversions_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)
    no_index_path = utils.mkstemp(dir=temp_dir)

    # Index will be used later for fast splitting
    # Will contain a tuple of (file position, number of lines) for every read
    # that has at least one mismatch.
    index = []
    no_index = []

    paired = {}  # Dictionary to store paired read information

    required_tags = []
    if not velocity:
        required_tags.append('GX')
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    # Can't use enumerate because n needs to be access outside from
    # the for loop.
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam, \
        open(conversions_path, 'w') as conversions_out, \
        open(no_conversions_path, 'w') as no_conversions_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = conversions_out.tell()
            read_id = read.query_name

            # Update every some interval. Updating every loop is very slow.
            n += 1
            if n % update_every == 0:
                lock.acquire()
                counter.value += update_every
                lock.release()

                # Update gene order (remove genes we no longer have to consider)
                # We take advantage of the fact that gene_order contains genes sorted
                # by coordinate.
                # When there are unmatched paired reads, we need to select the minimum start position
                # among those reads.
                if velocity:
                    start = read.reference_start if not paired else min(values[4] for values in paired.values())
                    to_remove = 0
                    for gene in gene_order:
                        if gene_infos[gene]['segment'].end <= start:
                            to_remove += 1
                        else:
                            break
                    if to_remove > 0:
                        gene_order = gene_order[to_remove:]

            # Skip reads that do not contain these tags
            if any(not read.has_tag(tag) for tag in required_tags) or read.is_duplicate or read.is_unmapped:
                continue

            barcode = read.get_tag(barcode_tag) if barcode_tag else 'NA'
            # If barcodes are provided, skip any barcodes not in the list.
            if barcodes and barcode not in barcodes:
                continue

            umi = read.get_tag(umi_tag) if umi_tag else 'NA'
            gx_assigned = read.has_tag('GX')
            gx = read.get_tag('GX') if gx_assigned else ''
            sequence = read.seq.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()
            reference_positions = read.get_reference_positions()

            # Count number of nucleotides in the region this read mapped to
            counts = Counter(reference)

            # Used exclusively for paired reads. So that we don't count conversions
            # in overlapping regions twice.
            # key: genome_i, value: (quality, original, converted)
            conversions = {}

            # If we are dealing with paired-end reads, we have to check for
            # overlaps.
            if read.is_paired:
                # Use HI tag, which is the multiple alignment index
                key = (read_id, read.get_tag('HI'))
                if key not in paired:
                    paired[key] = (
                        set(reference_positions),
                        reference,
                        sequence,
                        qualities,
                        read.reference_start,
                        read.get_aligned_pairs(matches_only=True, with_seq=True),
                    )
                    continue

                (
                    mate_reference_positions,  # Note that this is a set, not a list
                    mate_reference,
                    mate_sequence,
                    mate_qualities,
                    mate_reference_start,
                    aligned_pairs,
                ) = paired[key]
                del paired[key]

                # Update counts. union(A, B) = A + B - intersection(A, B)
                counts += Counter(mate_reference)
                overlap_pos = list(mate_reference_positions.intersection(reference_positions))
                overlap = Counter()
                for base, i in zip(reference, reference_positions):
                    if i in overlap_pos:
                        overlap[base] += 1
                counts -= overlap

                # Iterate through all aligned bases, but don't immediately write
                # to output. Instead, save the information in the conversions
                # dictionary because we need to combine this information with
                # the mate.
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

                # Update positions to be the union of the aligned position of the
                # mates. This is considered as one "read" to assign velocity.
                if velocity:
                    reference_positions = sorted(mate_reference_positions.union(reference_positions))

            # Assign velocity
            assignment = 'unassigned'
            if velocity:
                read_strand = None
                if strand == 'forward':
                    read_strand = '-' if read.is_reverse else '+'
                elif strand == 'reverse':
                    read_strand = '+' if read.is_reverse else '-'
                gx, assignment = assign_velocity_type(reference_positions, gx=gx, read_strand=read_strand)
                if not assignment or not gx:  # discard read if assignment failed
                    continue

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

                    conversions_out.write(
                        f'{read_id},{barcode},{umi},{gx},{contig},{genome_i},'
                        f'{genome_base},{read_base},{quality},'
                        f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]},'
                        f'{assignment},{gx_assigned}\n'
                    )

            # Flush cached conversions
            for genome_i, (quality, genome_base, read_base) in conversions.items():
                n_lines += 1
                conversions_out.write(
                    f'{read_id},{barcode},{umi},{gx},{contig},{genome_i},'
                    f'{genome_base},{read_base},{quality},'
                    f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]},'
                    f'{assignment},{gx_assigned}\n'
                )

            # Add to index if lines were written
            if n_lines > 0:
                index.append((pos, n_lines))
            else:
                no_index.append(no_conversions_out.tell())
                # Otherwise, this read did not contain any conversions.
                no_conversions_out.write(
                    f'{read_id},{barcode},{umi},{gx},'
                    f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]},'
                    f'{assignment},{gx_assigned}\n'
                )

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    # Save index
    index_path = utils.write_pickle(index, index_path, protocol=4)
    no_index_path = utils.write_pickle(no_index, no_index_path, protocol=4)

    return conversions_path, index_path, no_conversions_path, no_index_path


def check_bam_tags_exist(bam_path, tags, n_reads=100000, n_threads=8):
    """Utility function to check if BAM tags exists in a BAM within the first
    `n_reads` reads.

    :param bam_path: path to BAM
    :type bam_path: str
    :param tags: tags to check for
    :type tags: list
    :param n_reads: number of reads to consider, defaults to `100000`
    :type n_reads: int, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional

    :return: `True` if all tags were found, `False` otherwise
    :rtype: bool
    """
    tags_found = {tag: False for tag in tags}
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for i, read in enumerate(bam.fetch()):
            for tag in tags:
                if read.has_tag(tag):
                    tags_found[tag] = True

            if all(tags_found.values()):
                return True
            if i + 1 >= n_reads:
                break
    return False


def parse_all_reads(
    bam_path,
    conversions_path,
    index_path,
    no_conversions_path,
    no_index_path,
    gene_infos=None,
    transcript_infos=None,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    barcodes=None,
    n_threads=8,
    temp_dir=None,
    nasc=False,
    velocity=True
):
    """Parse all reads in a BAM and extract conversion, content and alignment
    information as CSVs.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions_path: path to output information about reads that have conversions
    :type conversions_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param no_conversions_path: path to output information about reads that do not have any conversions
    :type no_conversions_path: str
    :param no_index_path: path to no conversions index
    :type no_index_path: str
    :param gene_infos: dictionary containing gene information, as returned by
                       `preprocessing.gtf.parse_gtf`, defaults to `None`
    :type gene_infos: dictionary, optional
    :param transcript_infos: dictionary containing transcript information,
                             as returned by `preprocessing.gtf.parse_gtf`,
                             defaults to `None`
    :type transcript_infos: dictionary, optional
    :param strand: strandedness of the sequencing protocol, defaults to `forward`,
                   may be one of the following: `forward`, `reverse`, `None` (unstranded)
    :type strand: str, optional
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param barcodes: list of barcodes to be considered. All barcodes are considered
                     if not provided, defaults to `None`
    :type barcodes: list, optional
    :param n_threads: number of threads, defaults to `8`
    :type n_threads: int, optional
    :param temp_dir: path to temporary directory, defaults to `None`
    :type temp_dir: str, optional
    :param nasc: flag to change behavior to match NASC-seq pipeline, defaults to `False`
    :type nasc: bool, optional
    :param velocity: whether or not to assign a velocity type to each read,
                     defaults to `True`
    :type velocity: bool, optional

    :return: (path to conversions, path to conversions index,
              path to no conversions, path to no conversions index)
    :rtype: (str, str, str, str)
    """
    logger.debug('Checking if BAM has required tags')
    if not check_bam_tags_exist(bam_path, config.BAM_REQUIRED_TAGS, config.BAM_PEEK_READS, n_threads=n_threads):
        raise Exception(
            f'BAM does not contain required tags: {", ".join(config.BAM_REQUIRED_TAGS)}. '
            f'Please rerun alignment with at least `--outSAMattributes {" ".join(config.BAM_REQUIRED_TAGS)}`.'
        )

    logger.debug(f'Extracting contigs from BAM {bam_path}')
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            logger.debug(f'{index.total } reads mapped to contig `{index.contig}`')
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    if not velocity:
        args = [(contig, None, None) for contig in contigs]
    else:
        args = [(
            contig, {gene: info
                     for gene, info in gene_infos.items()
                     if info['chr'] == contig},
            {transcript: info
             for transcript, info in transcript_infos.items()
             if info['chr'] == contig}
        )
                for contig in contigs]

    # Initialize and run pool
    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            parse_read_contig,
            bam_path,
            counter,
            lock,
            strand=strand,
            umi_tag=umi_tag,
            barcode_tag=barcode_tag,
            barcodes=barcodes,
            temp_dir=tempfile.mkdtemp(dir=temp_dir),
            nasc=nasc,
            velocity=velocity,
        ), args
    )
    pool.close()

    # Display progres bar
    logger.debug(f'Processing contigs {contigs} from BAM')
    utils.display_progress_with_counter(counter, n_reads, async_result)
    pool.join()

    # Combine csvs
    logger.debug(f'Writing conversions to {conversions_path} and {no_conversions_path}')
    index = []
    no_index = []
    with open(conversions_path, 'wb') as conversions_out, \
        open(no_conversions_path, 'wb') as no_conversions_out:
        conversions_out.write(
            b'read_id,barcode,umi,GX,contig,genome_i,original,converted,quality,A,C,G,T,velocity,transcriptome\n'
        )
        no_conversions_out.write(b'read_id,barcode,umi,GX,A,C,G,T,velocity,transcriptome\n')

        pos = conversions_out.tell()
        no_pos = no_conversions_out.tell()

        for conversions_part_path, index_part_path, no_conversions_part_path, no_index_part_path in async_result.get():
            with open(conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, conversions_out)
            with open(no_conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, no_conversions_out)

            # Parse index
            index_part = utils.read_pickle(index_part_path)
            for p, n in index_part:
                index.append((pos + p, n))
            pos += os.path.getsize(conversions_part_path)

            no_index_part = utils.read_pickle(no_index_part_path)
            for p in no_index_part:
                no_index.append(no_pos + p)
            no_pos += os.path.getsize(no_conversions_part_path)
    logger.debug(f'Writing conversions index to {index_path} and {no_index_path}')
    index_path = utils.write_pickle(index, index_path, protocol=4)
    no_index_path = utils.write_pickle(no_index, no_index_path, protocol=4)

    return conversions_path, index_path, no_conversions_path, no_index_path
