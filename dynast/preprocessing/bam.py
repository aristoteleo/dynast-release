import multiprocessing
import os
import shutil
import tempfile
from collections import Counter
from functools import partial
from typing import List, Optional, Set, Tuple, Union

import ngs_tools as ngs
import numpy as np
import pandas as pd
import pysam
from typing_extensions import Literal

from .. import config, utils
from ..logging import logger

CONVERSION_CSV_COLUMNS = ['read_id', 'index', 'contig', 'genome_i', 'conversion', 'quality']
ALIGNMENT_COLUMNS = [
    'read_id', 'index', 'barcode', 'umi', 'GX', 'A', 'C', 'G', 'T', 'velocity', 'transcriptome', 'score'
]


def read_alignments(alignments_path: str, *args, **kwargs) -> pd.DataFrame:
    """Read alignments CSV as a pandas DataFrame.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    Args:
        alignments_path: path to alignments CSV

    Returns:
        Conversions dataframe
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


def read_conversions(conversions_path: str, *args, **kwargs) -> pd.DataFrame:
    """Read conversions CSV as a pandas DataFrame.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    Args:
        conversions_path: Path to conversions CSV

    Returns:
        Conversions dataframe
    """
    df = pd.read_csv(
        conversions_path,
        dtype={
            'read_id': 'string',
            'index': np.uint8,
            'contig': 'category',
            'genome_i': np.uint32,
            'conversion': 'category',
            'quality': np.uint8
        },
        na_filter=False,
        *args,
        **kwargs
    )
    return df


def select_alignments(df_alignments: pd.DataFrame) -> Set[Tuple[str, str]]:
    """Select alignments among duplicates. This function performs preliminary
    deduplication and returns a list of tuples (read_id, alignment index) to
    use for coverage calculation and SNP detection.

    Args:
        df_alignments: Alignments dataframe

    Returns:
        Set of (read_id, alignment index) that were selected
    """
    df_sorted = df_alignments.sort_values(['transcriptome', 'score'])
    umi = all(df_alignments['umi'] != 'NA')
    if umi:
        df_deduplicated = df_sorted.drop_duplicates(subset=['barcode', 'umi', 'GX'], keep='last')
    else:
        # TODO: This part can probably be removed because BAM parsing only considers
        # primary alignments now.

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
        counter: multiprocessing.Value,
        lock: multiprocessing.Lock,
        bam_path: str,
        contig: str,
        gene_infos: Optional[dict] = None,
        transcript_infos: Optional[dict] = None,
        strand: Literal['forward', 'reverse', 'unstranded'] = 'forward',
        umi_tag: Optional[str] = None,
        barcode_tag: Optional[str] = None,
        gene_tag: str = 'GX',
        barcodes: Optional[List[str]] = None,
        temp_dir: Optional[str] = None,
        update_every: int = 2000,
        nasc: bool = False,
        velocity: bool = True,
        strict_exon_overlap: bool = False,
) -> Tuple[str, str, str]:
    """Parse all reads mapped to a contig, outputing conversion
    information as temporary CSVs. This function is designed to be called as a
    separate process.

    Args:
        counter: Counter that keeps track of how many reads have been processed
        lock: Semaphore for the `counter` so that multiple processes do not
            modify it at the same time
        bam_path: Path to alignment BAM file
        contig: Only reads that map to this contig will be processed
        gene_infos: Dictionary containing gene information, as returned by
            `preprocessing.gtf.parse_gtf`, required if `velocity=True`
        transcript_infos: Dictionary containing transcript information,
            as returned by `preprocessing.gtf.parse_gtf`, required if `velocity=True`
        strand: Strandedness of the sequencing protocol, defaults to `forward`,
            may be one of the following: `forward`, `reverse`, `unstranded`
        umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
            `umi` column
        barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
            is output in the `barcode` column
        gene_tag: BAM tag that encodes gene assignment
        barcodes: List of barcodes to be considered. All barcodes are considered
            if not provided
        temp_dir: Path to temporary directory
        update_every: Update the counter every this many reads, defaults to `5000`
        nasc: Flag to change behavior to match NASC-seq pipeline, defaults to `False`
        velocity: Whether or not to assign a velocity type to each read
        strict_exon_overlap: Whether to use a stricter algorithm to assin reads
            as spliced

    Returns:
        (path to conversions, path to conversions index, path to alignments)
    """

    def skip_alignment(read, tags):
        return read.is_secondary or read.is_unmapped or read.is_duplicate or any(not read.has_tag(tag) for tag in tags)

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

            if any_exon_only and (not strict_exon_overlap or not any_intron_overlap):
                assigned_gene, assigned_type = gene, 'spliced'
            elif any_intron_overlap and (not strict_exon_overlap or not any_exon_only):
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
            if n == update_every:
                with lock:
                    counter.value += update_every
                n = 0

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

            # Skip reads
            if skip_alignment(read, required_tags):
                continue

            barcode = read.get_tag(barcode_tag) if barcode_tag else 'NA'
            if barcode == '-' or (barcodes and barcode not in barcodes):
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
                mate_sequence = mate.query_sequence.upper()
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
            sequence = read.query_sequence.upper()
            qualities = read.query_qualities
            reference = read.get_reference_sequence().upper()
            counts += Counter(reference)

            # Assign velocity
            assignment = 'unassigned'
            if velocity:
                read_strand = None
                if read.is_paired:
                    if read.is_read1:  # R1 is mapped after R2
                        if strand == 'forward':
                            read_strand = '+' if read.is_reverse else '-'
                        elif strand == 'reverse':
                            read_strand = '-' if read.is_reverse else '+'
                    else:  # R1 is mapped before R2
                        if strand == 'forward':
                            read_strand = '-' if read.is_reverse else '+'
                        elif strand == 'reverse':
                            read_strand = '+' if read.is_reverse else '-'
                elif strand == 'forward':
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
                        f'{genome_base}{read_base},{quality}\n'
                    )

            # Flush cached conversions
            for genome_i, (quality, genome_base, read_base) in conversions.items():
                n_lines += 1
                conversions_out.write(
                    f'{read_id},{alignment_index},{contig},{genome_i},'
                    f'{genome_base}{read_base},{quality}\n'
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

    with lock:
        counter.value += n

    # Save index
    index_path = utils.write_pickle(index, index_path, protocol=4)

    if gene_infos:
        del gene_infos
    if transcript_infos:
        del transcript_infos
    if barcodes:
        del barcodes

    return conversions_path, index_path, alignments_path


def get_tags_from_bam(bam_path: str, n_reads: int = 100000, n_threads: int = 8) -> Set[str]:
    """Utility function to retrieve all read tags present in a BAM.

    Args:
        bam_path: Path to BAM
        n_reads: Number of reads to consider
        n_threads: Number of threads

    Returns:
        Set of all tags found
    """
    tags = []
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads) as bam:
        for i, read in enumerate(bam.fetch(until_eof=True)):
            for (tag, _) in read.get_tags():
                tags.append(tag)
            if i + 1 >= n_reads:
                break
    return set(tags)


def check_bam_tags_exist(bam_path: str,
                         tags: List[str],
                         n_reads: int = 100000,
                         n_threads: int = 8) -> Tuple[bool, List[str]]:
    """Utility function to check if BAM tags exists in a BAM within the first
    `n_reads` reads.

    Args:
        bam_path: Path to BAM
        tags: Tags to check for
        n_reads: Number of reads to consider
        n_threads: Number of threads

    Returns:
        (whether all tags were found, list of not found tags)
    """
    tags_found = {tag: False for tag in tags}
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads) as bam:
        for i, read in enumerate(bam.fetch()):
            for tag in tags:
                if read.has_tag(tag):
                    tags_found[tag] = True

            if all(tags_found.values()):
                return True, []
            if i + 1 >= n_reads:
                break
    return False, [tag for tag in tags_found if not tags_found[tag]]


def check_bam_is_paired(bam_path: str, n_reads: int = 100000, n_threads: int = 8) -> bool:
    """Utility function to check if BAM has paired reads.

        bam_path: Path to BAM
        n_reads: Number of reads to consider
        n_threads: Number of threads

    Returns:
        Whether paired reads were detected
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads) as bam:
        for i, read in enumerate(bam.fetch()):
            if read.is_paired:
                return True

            if i + 1 >= n_reads:
                break
    return False


def check_bam_contains_secondary(bam_path: str, n_reads: int = 100000, n_threads: int = 8) -> bool:
    """Check whether BAM contains secondary alignments.

        bam_path: Path to BAM
        n_reads: Number of reads to consider
        n_threads: Number of threads

    Returns:
        Whether secondary alignments were detected
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads) as bam:
        for i, read in enumerate(bam.fetch()):
            if read.is_secondary:
                return True

            if i + 1 >= n_reads:
                break
    return False


def check_bam_contains_unmapped(bam_path: str) -> bool:
    """Check whether BAM contains unmapped reads.

        bam_path: Path to BAM

    Returns:
        Whether unmapped reads were detected
    """
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            if index.unmapped > 0:
                return True
    return False


def check_bam_contains_duplicate(bam_path, n_reads=100000, n_threads=8) -> bool:
    """Check whether BAM contains duplicates.

        bam_path: Path to BAM
        n_reads: Number of reads to consider
        n_threads: Number of threads

    Returns:
        Whether duplicates were detected
    """
    with pysam.AlignmentFile(bam_path, 'rb', threads=n_threads) as bam:
        for i, read in enumerate(bam.fetch()):
            if read.is_duplicate:
                return True

            if i + 1 >= n_reads:
                break
    return False


def sort_and_index_bam(bam_path: str, out_path: str, n_threads: int = 8, temp_dir: Optional[str] = None) -> str:
    """Sort and index BAM.

    If the BAM is already sorted, the sorting step is skipped.

    Args:
        bam_path: Path to alignment BAM file to sort
        out_path: Path to output sorted BAM
        n_threads: Number of threads
        temp_dir: Path to temporary directory

    Returns:
        Path to sorted and indexed BAM
    """
    bam_sorted = False
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        if f.header.get('HD', {}).get('SO') == 'coordinate':
            bam_sorted = True

    if not bam_sorted:
        logger.info(f'Sorting {bam_path} with samtools to {out_path}')
        args = ['-o', out_path, '-@', str(n_threads)]
        if temp_dir:
            args.extend(['-T', temp_dir])
        pysam.sort(bam_path, *args)
        bam_path = out_path

    # Check if BAM index exists and create one if it doesn't.
    bai_path = f'{bam_path}.bai'
    if not utils.all_exists(bai_path) or not bam_sorted:
        logger.info(f'Indexing {bam_path} with samtools to {bai_path}')
        pysam.index(bam_path, bai_path, '-@', str(n_threads))

    return bam_path


def split_bam(bam_path: str, n: int, n_threads: int = 8, temp_dir: Optional[str] = None) -> List[Tuple[str, int]]:
    """Split BAM into n parts.

    Args:
        bam_path: Path to alignment BAM file
        n: Number of splits
        n_threads: Number of threads
        temp_dir: Path to temporary directory

    Returns:
        List of tuples containing (split BAM path, number of reads)
    """
    # Filter for properly paired reads if BAM contains paired reads
    if check_bam_is_paired(bam_path, config.BAM_PEEK_READS, n_threads=n_threads):
        logger.debug('BAM contains paired reads. Filtering out non-properly paired reads')
        bam_path = ngs.bam.filter_bam(
            bam_path,
            lambda read: read.is_proper_pair if read.is_paired else True,
            utils.mkstemp(temp_dir, delete=True),
            n_threads=n_threads,
            show_progress=False,
        )
        pysam.index(bam_path, f'{bam_path}.bai', '-@', str(n_threads))

    logger.debug(f'Splitting BAM into {n} parts')
    splits = list(
        ngs.bam.split_bam(
            bam_path,
            utils.mkstemp(temp_dir, delete=True),
            n=n,
            n_threads=n_threads,
            show_progress=True,
        ).values()
    )
    for path, _ in splits:
        pysam.index(path, f'{path}.bai', '-@', str(n_threads))
    return splits


def parse_all_reads(
        bam_path: str,
        conversions_path: str,
        alignments_path: str,
        index_path: str,
        gene_infos: dict,
        transcript_infos: dict,
        strand: Literal['forward', 'reverse', 'unstranded'] = 'forward',
        umi_tag: Optional[str] = None,
        barcode_tag: Optional[str] = None,
        gene_tag: str = 'GX',
        barcodes: Optional[List[str]] = None,
        n_threads: int = 8,
        temp_dir: Optional[str] = None,
        nasc: bool = False,
        control: bool = False,
        velocity: bool = True,
        strict_exon_overlap: bool = False,
        return_splits: bool = False,
) -> Union[Tuple[str, str, str], Tuple[str, str, str, List[Tuple[str, int]]]]:
    """Parse all reads in a BAM and extract conversion, content and alignment
    information as CSVs.

        bam_path: Path to alignment BAM file
        conversions_path: Path to output information about reads that have conversions
        alignments_path: Path to alignments information about reads
        index_path: Path to conversions index
        no_index_path: Path to no conversions index
        gene_infos: Dictionary containing gene information, as returned by
            `ngs.gtf.genes_and_transcripts_from_gtf`
        transcript_infos: Dictionary containing transcript information,
            as returned by `ngs.gtf.genes_and_transcripts_from_gtf`
        strand: Strandedness of the sequencing protocol, defaults to `forward`,
            may be one of the following: `forward`, `reverse`, `unstranded`
        umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
            `umi` column
        barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
            is output in the `barcode` column
        gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
        barcodes: List of barcodes to be considered. All barcodes are considered
            if not provided
        n_threads: Number of threads
        temp_dir: Path to temporary directory
        nasc: Flag to change behavior to match NASC-seq pipeline
        velocity: Whether or not to assign a velocity type to each read
        strict_exon_overlap: Whether to use a stricter algorithm to assign reads as spliced
        return_splits: return BAM splits for later reuse

    Returns:
        (path to conversions, path to alignments, path to conversions index)
            If `return_splits` is True, then there is an additional return value, which
            is a list of tuples containing split BAM paths and number of reads
            in each BAM.
    """
    contig_genes = {}
    contig_transcripts = {}
    for gene_id, gene_info in gene_infos.items():
        contig_genes.setdefault(gene_info['chromosome'], []).append(gene_id)
        contig_transcripts.setdefault(gene_info['chromosome'], []).extend(gene_info['transcripts'])

    splits = split_bam(
        bam_path, min(n_threads * 8,
                      utils.get_file_descriptor_limit() - 10), n_threads=n_threads, temp_dir=temp_dir
    ) if n_threads > 1 else [(bam_path, ngs.bam.count_bam(bam_path))]

    # Add argument per contig
    args = []
    for path, _ in splits:
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
            strict_exon_overlap=strict_exon_overlap,
        ), args
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, sum(value[1] for value in splits), async_result)
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

    return (conversions_path, alignments_path,
            index_path) if not return_splits else (conversions_path, alignments_path, index_path, splits)
