import multiprocessing
import re
from collections import Counter
from functools import partial
from typing import Dict, List, Optional, Set, Tuple

import pysam

from .. import utils
from ..logging import logger

COVERAGE_PARSER = re.compile(
    r'''^
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<coverage>[^,]*)\n
    $''', re.VERBOSE
)


def read_coverage(coverage_path: str) -> Dict[str, Dict[int, int]]:
    """Read coverage CSV as a dictionary.

    Args:
        coverage_path: Path to coverage CSV

    Returns:
        Coverage as a nested dictionary
    """
    coverage = {}
    with open(coverage_path, 'r') as f:
        for i, line in enumerate(f):
            groups = COVERAGE_PARSER.match(line).groupdict()
            try:
                contig, genome_i, count = groups['contig'], int(groups['genome_i']), int(groups['coverage'])
                coverage.setdefault(contig, {})[genome_i] = count
            except ValueError as e:
                if i == 0:
                    continue
                raise e
    return coverage


def calculate_coverage_contig(
        counter: multiprocessing.Value,
        lock: multiprocessing.Lock,
        bam_path: str,
        contig: str,
        indices: List[Tuple[int, int, int]],
        alignments: Set[Tuple[str, int]] = None,
        umi_tag: Optional[str] = None,
        barcode_tag: Optional[str] = None,
        gene_tag: str = 'GX',
        barcodes: Optional[List[str]] = None,
        temp_dir: Optional[str] = None,
        update_every: int = 50000,
        velocity: bool = True
) -> str:
    """Calculate converage for a specific contig. This function is designed to
    be called as a separate process.

    Args:
        counter: Counter that keeps track of how many reads have been processed
        lock: Semaphore for the `counter` so that multiple processes do not
            modify it at the same time
        bam_path: Path to alignment BAM file
        contig: Only reads that map to this contig will be processed
        indices: Genomic positions to consider
        alignments: Set of (read_id, alignment_index) tuples to process. All
            alignments are processed if this option is not provided.
        umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
            `umi` column
        barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
            is output in the `barcode` column
        gene_tag: BAM tag that encodes gene assignment, defaults to `GX`
        barcodes: List of barcodes to be considered. All barcodes are considered if not provided
        temp_dir: Path to temporary directory
        update_every: Update the counter every this many reads
        velocity: Whether or not velocities were assigned

    Returns:
        Path to coverage CSV
    """

    def skip_alignment(read, tags):
        return read.is_secondary or read.is_unmapped or read.is_duplicate or any(not read.has_tag(tag) for tag in tags)

    coverage_path = utils.mkstemp(dir=temp_dir)
    paired = {}
    coverage = Counter()

    required_tags = []
    if not velocity:
        required_tags.append(gene_tag)
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    n = 0
    with pysam.AlignmentFile(bam_path, 'rb', threads=2) as bam, \
        open(coverage_path, 'w') as coverage_out:
        for read in bam.fetch(contig):
            n += 1
            if n == update_every:
                with lock:
                    counter.value += update_every
                n = 0

            if len(coverage) > 1000000:
                # NOTE: dictionary keys are sorted by insertion order
                for genome_i in list(coverage.keys()):
                    if genome_i < read.reference_start:
                        coverage_out.write(f'{contig},{genome_i},{coverage[genome_i]}\n')
                        del coverage[genome_i]

                    else:
                        break

            if skip_alignment(read, required_tags):
                continue

            barcode = read.get_tag(barcode_tag) if barcode_tag else 'NA'
            if barcode == '-' or (barcodes and barcode not in barcodes):
                continue

            read_id = read.query_name
            alignment_index = read.get_tag('HI')
            key = (read_id, alignment_index)
            if alignments and key not in alignments:
                continue

            reference_positions = read.get_reference_positions()

            # For paired end reads, don't count overlap twice
            if read.is_paired:
                if key not in paired:
                    paired[key] = reference_positions
                    continue

                reference_positions = list(set(reference_positions + paired[key]))
                del paired[key]

            for genome_i in reference_positions:
                if genome_i in indices:
                    coverage[genome_i] += 1

        for genome_i, cover in coverage.items():
            coverage_out.write(f'{contig},{genome_i},{cover}\n')

    with lock:
        counter.value += n
    del coverage
    del indices
    if alignments:
        del alignments
    if barcodes:
        del barcodes

    return coverage_path


def calculate_coverage(
        bam_path: str,
        conversions: Dict[str, Set[int]],
        coverage_path: str,
        alignments: Optional[List[Tuple[str, int]]] = None,
        umi_tag: Optional[str] = None,
        barcode_tag: Optional[str] = None,
        gene_tag: str = 'GX',
        barcodes: Optional[List[str]] = None,
        temp_dir: Optional[str] = None,
        velocity: bool = True,
) -> str:
    """Calculate coverage of each genomic position per barcode.

    Args:
        bam_path: Path to alignment BAM file
        conversions: Dictionary of contigs as keys and sets of genomic positions
            as values that indicates positions where conversions were observed
        coverage_path: Path to write coverage CSV
        alignments: Set of (read_id, alignment_index) tuples to process. All
            alignments are processed if this option is not provided.
        umi_tag: BAM tag that encodes UMI, if not provided, `NA` is output in the
            `umi` column
        barcode_tag: BAM tag that encodes cell barcode, if not provided, `NA`
            is output in the `barcode` column
        gene_tag: BAM tag that encodes gene assignment
        barcodes: List of barcodes to be considered. All barcodes are considered if not provided
        temp_dir: Path to temporary directory
        velocity: Whether or not velocities were assigned

    Returns:
        Path to coverage CSV
    """
    logger.debug(f'Extracting contigs from BAM {bam_path}')
    contigs = []
    n_reads = 0
    with pysam.AlignmentFile(bam_path, 'rb') as bam:
        for index in bam.get_index_statistics():
            n_reads += index.total
            if index.total > 0:
                contigs.append(index.contig)

    # Turns out coverage calculation per read is fast enough that we end
    # up being IO-bound if we use multiprocessing. So, we just process
    # sequentially (per contig) here instead.
    # We create a new pool because we still want to be able to keep track of
    # the current progress.
    pool, counter, lock = utils.make_pool_with_counter(1)
    async_result = pool.starmap_async(
        partial(
            calculate_coverage_contig,
            counter,
            lock,
            bam_path,
            alignments=alignments,
            umi_tag=umi_tag,
            barcode_tag=barcode_tag,
            gene_tag=gene_tag,
            barcodes=barcodes,
            velocity=velocity
        ), [(contig, conversions.get(contig, set())) for contig in contigs]
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, n_reads, async_result)
    pool.join()

    logger.debug('Combining results')
    coverage = {}
    for coverage_part_path in async_result.get():
        cov = read_coverage(coverage_part_path)
        coverage = utils.merge_dictionaries(coverage, cov)

    logger.debug(f'Writing coverage to {coverage_path}')
    with open(coverage_path, 'w') as coverage_out:
        coverage_out.write('contig,genome_i,coverage\n')
        for (contig, genome_i), cover in utils.flatten_dictionary(coverage):
            coverage_out.write(f'{contig},{genome_i},{cover}\n')

    return coverage_path
