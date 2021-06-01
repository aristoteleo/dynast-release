import os
import shutil
import tempfile
from collections import Counter
from functools import partial

import ngs_tools as ngs
import numpy as np
import pandas as pd
import pysam

from .. import config, utils
from ..logging import logger

CONVERSION_CSV_COLUMNS = ['read_id', 'index', 'contig', 'genome_i', 'original', 'converted', 'quality']
ALIGNMENT_COLUMNS = [
    'read_id', 'index', 'barcode', 'umi', 'GX', 'A', 'C', 'G', 'T', 'velocity', 'transcriptome', 'score'
]


def read_alignments(alignments_path, *args, **kwargs):
    """Read alignments CSV as a pandas DataFrame.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    :param alignments_path: path to alignments CSV
    :type alignments_path: str

    :return: conversions dataframe
    :rtype: pandas.DataFrame
    """
    df = pd.read_csv(
        alignments_path,
        dtype={
            'read_id': 'string',
            'index': np.uint8,
            'barcode': 'category',
            'umi': 'category',
            'GX': 'category',
            'A': np.uint8,
            'C': np.uint8,
            'G': np.uint8,
            'T': np.uint8,
            'velocity': 'category',
            'transcriptome': bool,
            'score': np.uint16
        },
        na_filter=False,
        *args,
        **kwargs
    )
    return df


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
            'index': np.uint8,
            'contig': 'category',
            'genome_i': np.uint32,
            'original': 'category',
            'converted': 'category',
        },
        na_filter=False,
        *args,
        **kwargs
    )
    return df


def select_alignments(df_alignments):
    """Select alignments among duplicates. This function performs preliminary
    deduplication and returns a list of tuples (read_id, alignment index) to
    use for coverage calculation and SNP detection.

    :param df_alignments: alignments dataframe
    :type df_alignments: pandas.DataFrame

    :return: set of (read_id, alignment index) that were selected
    :rtype: set
    """
    df_sorted = df_alignments.sort_values(['transcriptome', 'score'])
    umi = all(df_alignments['umi'] != 'NA')
    if umi:
        df_deduplicated = df_sorted.drop_duplicates(subset=['barcode', 'umi', 'GX'], keep='last')
    else:
        # Drop any multimappers where none map to the transcriptome and
        # are assigned multiple genes
        df_sorted['not_transcriptome'] = ~df_sorted['transcriptome']
        read_id_grouped = df_sorted.groupby('read_id', sort=False, observed=True)
        not_transcriptome = read_id_grouped['not_transcriptome'].all()
        not_transcriptome_read_ids = not_transcriptome.index[not_transcriptome]

        multigene = read_id_grouped['GX'].nunique() > 1
        multigene_read_ids = multigene.index[multigene]

        to_drop = list(not_transcriptome_read_ids.intersection(multigene_read_ids))
        df_deduplicated = df_sorted[~df_sorted['read_id'].isin(to_drop)].drop_duplicates('read_id', keep='last')

    return set(df_deduplicated[['read_id', 'index']].itertuples(index=False, name=None))


def parse_read_contig(
    counter,
    lock,
    bam_path,
    contig,
    gene_infos=None,
    transcript_infos=None,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    barcodes=None,
    temp_dir=None,
    update_every=2000,
    nasc=False,
    velocity=True
):
    """Parse all reads mapped to a contig, outputing conversion
    information as temporary CSVs. This function is designed to be called as a
    separate process.

    :param counter: counter that keeps track of how many reads have been processed
    :type counter: multiprocessing.Value
    :param lock: semaphore for the `counter` so that multiple processes do not
                 modify it at the same time
    :type lock: multiprocessing.Lock
    :param bam_path: path to alignment BAM file
    :type bam_path: str
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
    :param gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
    :type gene_tag: str, optional
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

    :return: (path to conversions, path to conversions index, path to alignments)
    :rtype: (str, str, str)
    """
    # Sort genes by segment positions
    if velocity:
        gene_order = sorted(gene_infos.keys(), key=lambda gene: tuple(gene_infos[gene]['segment']))

    def assign_velocity_type(positions, gx=None, read_strand=None):
        alignment = ngs.gtf.SegmentCollection.from_positions(positions)

        # If gx is provided, narrow down search to this single gene
        assigned_gene = gx
        assigned_type = None
        for gene in gene_order if not gx else [gx]:
            # Skip or break loop conditions, used when gx is not provided
            if not gx:
                # On different strand
                if read_strand and read_strand != gene_infos[gene]['strand']:
                    continue
                # Out of range
                gene_segment = gene_infos[gene]['segment']
                if not (alignment.start >= gene_segment.start and alignment.end <= gene_segment.end):
                    continue
                elif alignment.end <= gene_segment.start:
                    break

            any_exon_only = False
            any_exon_overlap = False
            any_intron_overlap = False
            for transcript_id in gene_infos[gene]['transcripts']:
                transcript = transcript_infos[transcript_id]

                if not any_exon_overlap:
                    any_exon_overlap |= alignment.is_overlapping(transcript['exons'])
                if any_exon_overlap and not any_exon_only:
                    any_exon_only |= alignment.is_subset(transcript['exons'])
                if not any_intron_overlap:
                    any_intron_overlap |= alignment.is_overlapping(transcript['introns'])

            # Additional conditions for when gx is not provided
            if not gx:
                # Multi-gene
                if assigned_gene and (any_exon_overlap or any_intron_overlap):
                    return None, None
                # No overlaps
                if not any_exon_overlap and not any_intron_overlap:
                    continue

            if any_exon_only and not any_intron_overlap:
                assigned_gene, assigned_type = gene, 'spliced'
            elif not any_exon_only and any_intron_overlap:
                assigned_gene, assigned_type = gene, 'unspliced'
            else:
                assigned_gene, assigned_type = gene, 'ambiguous'

        return assigned_gene, assigned_type

    conversions_path = utils.mkstemp(dir=temp_dir)
    alignments_path = utils.mkstemp(dir=temp_dir)
    index_path = utils.mkstemp(dir=temp_dir)

    # Index will be used later for fast splitting
    # Will contain a tuple of (file position, number of lines) for every read
    # that has at least one mismatch.
    index = []

    paired = {}  # Dictionary to store paired read information

    # Any reads that don't have all the required tags will be skipped.
    required_tags = []
    if not velocity:
        required_tags.append(gene_tag)
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    # Can't use enumerate because n needs to be access outside from
    # the for loop.
    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(conversions_path, 'w') as conversions_out, \
        open(alignments_path, 'w') as alignments_out:
        for read in bam.fetch(contig):
            n_lines = 0
            pos = conversions_out.tell()

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
                    # IMPORTANT: as of Python 3.7, dictionary order is guaranteed to be
                    # insertion order.
                    # https://docs.python.org/3.7/library/stdtypes.html#dict
                    start = read.reference_start if not paired else next(iter(paired.values())).reference_start
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

            read_id = read.query_name
            alignment_index = read.get_tag('HI')
            alignment_score = read.get_tag('AS')
            reference_positions = read.get_reference_positions()
            counts = Counter()

            # Used exclusively for paired reads. So that we don't count conversions
            # in overlapping regions twice.
            # key: genome_i, value: (quality, original, converted)
            conversions = {}

            # If we are dealing with paired-end reads, we have to check for
            # overlaps.
            if read.is_paired:
                # Use HI tag, which is the multiple alignment index
                key = (read_id, alignment_index)
                if key not in paired:
                    paired[key] = read
                    continue

                mate = paired[key]
                mate_sequence = mate.seq.upper()
                mate_qualities = mate.query_qualities
                mate_reference = mate.get_reference_sequence().upper()
                mate_reference_positions = mate.get_reference_positions()
                del paired[key]

                # Update counts. union(A, B) = A + B - intersection(A, B)
                counts += Counter(mate_reference)
                overlap_pos = list(set(mate_reference_positions).intersection(reference_positions))
                overlap = Counter()

                # Iterate through all aligned bases, but don't immediately write
                # to output. Instead, save the information in the conversions
                # dictionary because we need to combine this information with
                # the mate.
                for read_i, genome_i, _genome_base in mate.get_aligned_pairs(matches_only=True, with_seq=True):
                    read_base = mate_sequence[read_i]
                    genome_base = _genome_base.upper()
                    if 'N' in (genome_base, read_base):
                        continue

                    if genome_i in overlap_pos:
                        overlap[genome_base] += 1

                    if genome_base != read_base and (not nasc or genome_i not in overlap_pos):
                        # Add this conversion to the conversions dictionary.
                        # We don't write this conversion immediately because
                        # we want to take the conversion with higher quality
                        # if conversions exist in same positions in overlapping
                        # segments.
                        conversions[genome_i] = (mate_qualities[read_i], genome_base, read_base)
                counts -= overlap
                del mate

                # Update positions to be the union of the aligned position of the
                # mates. This is considered as one "read" to assign velocity.
                if velocity:
                    reference_positions += mate_reference_positions

            umi = read.get_tag(umi_tag) if umi_tag else 'NA'
            gx_assigned = read.has_tag(gene_tag)
            gx = read.get_tag(gene_tag) if gx_assigned else ''
            sequence = read.seq.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()
            counts += Counter(reference)

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
                        # The mate's conversion has higher quality, so replace
                        # output with that of this mate.
                        if mate_quality > quality:
                            read_base = mate_read_base
                            quality = mate_quality
                    del conversions[genome_i]

                if genome_base != read_base:
                    n_lines += 1

                    conversions_out.write(
                        f'{read_id},{alignment_index},{contig},{genome_i},'
                        f'{genome_base},{read_base},{quality}\n'
                    )

            # Flush cached conversions
            for genome_i, (quality, genome_base, read_base) in conversions.items():
                n_lines += 1
                conversions_out.write(
                    f'{read_id},{alignment_index},{contig},{genome_i},'
                    f'{genome_base},{read_base},{quality}\n'
                )

            # Add to index if lines were written
            if n_lines > 0:
                index.append((pos, n_lines, alignments_out.tell()))
            # Write alignment info
            alignments_out.write(
                f'{read_id},{alignment_index},{barcode},{umi},{gx},'
                f'{counts["A"]},{counts["C"]},{counts["G"]},{counts["T"]},'
                f'{assignment},{gx_assigned},{alignment_score}\n'
            )

    lock.acquire()
    counter.value += n % update_every
    lock.release()

    # Save index
    index_path = utils.write_pickle(index, index_path, protocol=4)

    return conversions_path, index_path, alignments_path


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

    :return: (whether all tags were found, list of not found tags)
    :rtype: (bool, list)
    """
    tags_found = {tag: False for tag in tags}
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for i, read in enumerate(bam.fetch()):
            for tag in tags:
                if read.has_tag(tag):
                    tags_found[tag] = True

            if all(tags_found.values()):
                return True, []
            if i + 1 >= n_reads:
                break
    return False, [tag for tag in tags_found if not tags_found[tag]]


def parse_all_reads(
    bam_path,
    conversions_path,
    alignments_path,
    index_path,
    gene_infos,
    transcript_infos,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    barcodes=None,
    n_threads=8,
    temp_dir=None,
    nasc=False,
    control=False,
    velocity=True
):
    """Parse all reads in a BAM and extract conversion, content and alignment
    information as CSVs.

    :param bam_path: path to alignment BAM file
    :type bam_path: str
    :param conversions_path: path to output information about reads that have conversions
    :type conversions_path: str
    :param alignments_path: path to alignments information about reads
    :type alignments_path: str
    :param index_path: path to conversions index
    :type index_path: str
    :param no_index_path: path to no conversions index
    :type no_index_path: str
    :param gene_infos: dictionary containing gene information, as returned by
                       `ngs.gtf.genes_and_transcripts_from_gtf`
    :type gene_infos: dictionary
    :param transcript_infos: dictionary containing transcript information,
                             as returned by `ngs.gtf.genes_and_transcripts_from_gtf`
    :type transcript_infos: dictionary
    :param strand: strandedness of the sequencing protocol, defaults to `forward`,
                   may be one of the following: `forward`, `reverse`, `unstranded`
    :type strand: str, optional
    :param umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
                    `umi` column, defaults to `None`
    :type umi_tag: str, optional
    :param barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
                        is output in the `barcode` column, defaults to `None`
    :type barcode_tag: str, optional
    :param gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
    :type gene_tag: str, optional
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

    :return: (path to conversions, path to alignments, path to conversions index)
    :rtype: (str, str, str)
    """
    logger.debug('Checking if BAM has required tags')
    tags = config.BAM_REQUIRED_TAGS.copy()
    if barcode_tag:
        tags.append(barcode_tag)
    if umi_tag:
        tags.append(umi_tag)
    if gene_tag:
        tags.append(gene_tag)
    all_exists, not_found = check_bam_tags_exist(bam_path, tags, config.BAM_PEEK_READS, n_threads=n_threads)
    if not all_exists:
        raise Exception(
            f'First {config.BAM_PEEK_READS} reads in the BAM do not contain the following required tags: '
            f'{", ".join(not_found)}. '
        )

    contig_genes = {}
    contig_transcripts = {}
    for gene_id, gene_info in gene_infos.items():
        contig_genes.setdefault(gene_info['chromosome'], []).append(gene_id)
        contig_transcripts.setdefault(gene_info['chromosome'], []).extend(gene_info['transcripts'])

    if n_threads > 1:
        n_splits = min(n_threads * 8, utils.get_file_descriptor_limit() - 10)
        logger.debug(f'Splitting BAM into {n_splits} parts')
        splits = ngs.bam.split_bam(
            bam_path,
            utils.mkstemp(temp_dir, delete=True),
            n=n_splits,
            n_threads=n_threads,
            show_progress=False,
        )
        for path, _ in splits.values():
            pysam.index(path, f'{path}.bai', '-@', str(n_threads))

    else:
        splits = {'0': (bam_path, ngs.bam.count_bam(bam_path))}

    # Add argument per contig
    args = []
    for path, _ in splits.values():
        with pysam.AlignmentFile(path, 'rb') as f:
            for index in f.get_index_statistics():
                if index.total > 0:
                    args.append((
                        path, index.contig,
                        {gene_id: gene_infos[gene_id]
                         for gene_id in contig_genes.get(index.contig, {})} if velocity else None, {
                             transcript_id: transcript_infos[transcript_id]
                             for transcript_id in contig_transcripts.get(index.contig, {})
                         } if velocity else None
                    ))

    # Initialize and run pool
    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    async_result = pool.starmap_async(
        partial(
            parse_read_contig,
            counter,
            lock,
            strand=strand,
            umi_tag=umi_tag,
            barcode_tag=barcode_tag,
            gene_tag=gene_tag,
            barcodes=barcodes,
            temp_dir=tempfile.mkdtemp(dir=temp_dir),
            nasc=nasc,
            velocity=velocity,
        ), args
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, sum(value[1] for value in splits.values()), async_result)
    pool.join()

    # Combine csvs
    logger.debug(f'Writing alignments to {alignments_path} and conversions to {conversions_path}')
    index = []
    with open(conversions_path, 'wb') as conversions_out, \
        open(alignments_path, 'wb') as alignments_out:
        conversions_out.write(f'{",".join(CONVERSION_CSV_COLUMNS)}\n'.encode())
        alignments_out.write(f'{",".join(ALIGNMENT_COLUMNS)}\n'.encode())

        pos = conversions_out.tell()
        pos2 = alignments_out.tell()

        for conversions_part_path, index_part_path, alignments_part_path in async_result.get():
            with open(conversions_part_path, 'rb') as f:
                shutil.copyfileobj(f, conversions_out)
            with open(alignments_part_path, 'rb') as f:
                shutil.copyfileobj(f, alignments_out)

            # Parse index
            index_part = utils.read_pickle(index_part_path)
            for p, n, a in index_part:
                index.append((pos + p, n, pos2 + a))
            pos += os.path.getsize(conversions_part_path)
            pos2 += os.path.getsize(alignments_part_path)
    logger.debug(f'Writing conversions index to {index_path}')
    index_path = utils.write_pickle(index, index_path)

    return conversions_path, alignments_path, index_path
