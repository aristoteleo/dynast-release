import array
import multiprocessing
import queue
from hashlib import sha256

import pysam
import ngs_tools as ngs
import numpy as np
from tqdm import tqdm

from . import bam
from .. import utils
from ..logging import logger

BASES = ('A', 'C', 'G', 'T')
BASE_IDX = {base: i for i, base in enumerate(BASES)}


def call_consensus_from_reads(reads, header, quality=27):
    """Call a single consensus alignment given a list of aligned reads.

    Reads must map to the same contig. Results are undefined otherwise.
    Additionally, consensus bases are called only for positions that match
    to the reference (i.e. no insertions allowed).

    This function only sets the minimal amount of attributes such that the
    alignment is valid. These include:
    * read name -- SHA256 hash of the provided read names
    * read sequence and qualities
    * reference name and ID
    * reference start
    * mapping quality (MAPQ)
    * cigarstring
    * MD tag
    * NM tag
    * Not unmapped, paired, duplicate, qc fail, secondary, nor supplementary

    The caller is expected to further populate the alignment
    with additional tags, flags, and name.

    :param reads: List of reads to call a consensus sequence from
    :type reads: list
    :param header: header to use when creating the new pysam alignment
    :type header: pysam.AlignmentHeader
    :param quality: quality threshold, defaults to 27
    :type quality: int, optional

    :return: (New pysam alignment of the consensus sequence)
    :rtype: pysam.AlignedSegment
    """
    if len(set(read.reference_name for read in reads)) > 1:
        raise Exception("Can not call consensus from reads mapping to multiple contigs.")

    # Pysam coordinates are [start, end)
    reads = sorted(reads, key=lambda read: read.reference_start)
    left_pos = reads[0].reference_start
    right_pos = max(read.reference_end for read in reads)
    length = right_pos - left_pos

    # A consensus sequence is internally represented as a L x 4 matrix,
    # where L is the length of the sequence and the columns correspond to
    # each of the four bases. The values indicate the support of each base.
    # It's possible to switch these to sparse matrices if memory becomes an issue.
    sequence = np.zeros((length, len(BASES)), dtype=np.uint32)
    reference = np.full(length, -1, dtype=np.int8)  # -1 means unobserved

    for read in reads:
        read_sequence = read.query_sequence.upper()
        read_qualities = read.query_qualities
        for read_i, genome_i, _genome_base in read.get_aligned_pairs(matches_only=True, with_seq=True):
            read_base = read_sequence[read_i]
            genome_base = _genome_base.upper()
            if 'N' in (genome_base, read_base):
                continue

            i = genome_i - left_pos
            if reference[i] < 0:
                reference[i] = BASE_IDX[genome_base]
            sequence[i, BASE_IDX[read_base]] += read_qualities[read_i]

    # Determine consensus
    # Note that we ignore any insertions and deletions.
    consensus_length = (reference >= 0).sum()
    consensus = np.zeros(consensus_length, dtype=np.uint8)
    qualities = np.zeros(consensus_length, dtype=np.uint8)
    cigar = []
    last_cigar_op = None
    cigar_n = 0
    md = []
    md_n = 0
    md_zero = True
    nm = 0
    consensus_i = 0
    for i in range(length):
        ref = reference[i]
        # Region not present in read. MD tag only deals with aligned
        # regions, so nothing else needs to be done.
        cigar_op = 'N'
        if ref >= 0:
            seq = sequence[i]
            base = seq.argmax()
            base_q = seq[base]
            if base_q < quality:
                base = ref

            # We use the STAR convention of using M cigar operation to mean
            # both matches AND mismatches, ignoring the X cigar operation exists.
            cigar_op = 'M'

            if ref == base:
                md_n += 1
                md_zero = False
            else:
                if md_n > 0 or md_zero:
                    md.append(str(md_n))
                    md_n = 0
                md.append(BASES[ref])
                md_zero = True
                nm += 1

            consensus[consensus_i] = base
            qualities[consensus_i] = min(base_q, 42)  # Clip to maximum PHRED score
            consensus_i += 1

        if cigar_op == last_cigar_op:
            cigar_n += 1
        else:
            if last_cigar_op:
                cigar.append(f'{cigar_n}{last_cigar_op}')
            last_cigar_op = cigar_op
            cigar_n = 1

    md.append(str(md_n))  # MD tag always ends with a number
    cigar.append(f'{cigar_n}{last_cigar_op}')

    al = pysam.AlignedSegment(header)
    al.query_name = sha256(''.join(read.query_name for read in reads).encode('utf-8')).hexdigest()
    al.query_sequence = ''.join(BASES[i] for i in consensus)
    al.query_qualities = array.array('B', qualities)
    al.reference_name = reads[0].reference_name
    al.reference_id = reads[0].reference_id
    al.reference_start = left_pos
    al.mapping_quality = 255
    al.cigarstring = ''.join(cigar)
    al.set_tags([('MD', ''.join(md)), ('NM', nm)])

    # Make sure these are False
    al.is_unmapped = False
    al.is_paired = False
    al.is_duplicate = False
    al.is_qcfail = False
    al.is_secondary = False
    al.is_supplementary = False
    return al


def call_consensus_from_reads_process(reads, header, tags, strand=None, quality=27):
    header = pysam.AlignmentHeader.from_dict(header)
    reads = [pysam.AlignedSegment.fromstring(read, header) for read in reads]
    consensus = call_consensus_from_reads(reads, header, quality=quality)
    consensus.set_tags(tags)
    consensus.is_paired = False
    if strand == '-':
        consensus.is_reverse = True
    return consensus.to_string()


def consensus_worker(args_q, results_q, *args, **kwargs):
    while True:
        try:
            _args = args_q.get(timeout=1)  # None means we are done.
        except queue.Empty:
            continue
        if _args is None:
            return

        results_q.put(call_consensus_from_reads_process(*_args, *args, **kwargs))


def call_consensus(
    bam_path,
    out_path,
    gene_infos,
    strand='forward',
    umi_tag=None,
    barcode_tag=None,
    gene_tag='GX',
    quality=27,
    temp_dir=None,
    n_threads=8
):

    def skip_alignment(read, tags):
        return read.is_secondary or read.is_unmapped or any(not read.has_tag(tag) for tag in tags)

    def find_genes(contig, start, end, read_strand=None):
        genes = []
        for gene in contig_gene_order[contig]:
            if read_strand and read_strand != gene_infos[gene]['strand']:
                continue
            gene_segment = gene_infos[gene]['segment']
            if end <= gene_segment.start:
                break
            if start >= gene_segment.start and end <= gene_segment.end:
                genes.append(gene)
        return genes

    contig_gene_order = {}
    for gene_id, gene_info in gene_infos.items():
        contig_gene_order.setdefault(gene_info['chromosome'], []).append(gene_id)
    for contig in list(contig_gene_order.keys()):
        contig_gene_order[contig] = sorted(
            contig_gene_order[contig], key=lambda gene: tuple(gene_infos[gene]['segment'])
        )

    gx_barcode_umi_groups = {}
    paired = {}

    required_tags = []
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    # Start processes for consensus calling
    logger.debug(f'Spawning {n_threads} processes')
    manager = multiprocessing.Manager()
    args_q = manager.Queue(1000 * n_threads)
    results_q = manager.Queue()
    workers = [
        multiprocessing.Process(
            target=consensus_worker, args=(args_q, results_q), kwargs=dict(quality=quality), daemon=True
        ) for _ in range(n_threads)
    ]
    for worker in workers:
        worker.start()

    temp_out_path = utils.mkstemp(dir=temp_dir)
    with pysam.AlignmentFile(bam_path, 'rb') as f:
        # Get header dict and update sort order to unsorted.
        header_dict = f.header.to_dict()
        hd = header_dict.setdefault('HD', {'VN': '1.4', 'SO': 'unsorted'})
        hd['SO'] = 'unsorted'
        header = pysam.AlignmentHeader.from_dict(header_dict)
        with pysam.AlignmentFile(temp_out_path, 'wb', header=header) as out:
            for i, read in tqdm(enumerate(f.fetch()), total=ngs.bam.count_bam(bam_path), ascii=True, smoothing=0.01,
                                desc='Calling consensus'):
                if skip_alignment(read, required_tags):
                    out.write(read)
                    continue

                contig = read.reference_name
                barcode = read.get_tag(barcode_tag) if barcode_tag else None
                umi = read.get_tag(barcode_tag) if umi_tag else None

                read_id = read.query_name
                alignment_index = read.get_tag('HI')
                start = read.reference_start
                end = read.reference_end
                key = (read_id, alignment_index)
                mate = None
                if read.is_paired:
                    if key not in paired:
                        paired[key] = read
                        continue

                    mate = paired.pop(key)
                    # Use alignment start and end as UMI for paired reads without UMI
                    if not umi:
                        start = mate.reference_start
                        umi = (start, end)

                # Determine read strand
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

                # Find compatible genes
                gx_assigned = read.has_tag(gene_tag)
                genes = [read.get_tag(gene_tag)] if gx_assigned else find_genes(contig, start, end, read_strand)

                # If there isn't exactly one compatible gene, do nothing and
                # write to BAM.
                if len(genes) != 1:
                    out.write(read)
                    if read.is_paired:
                        out.write(mate)
                    continue

                # Add read to group
                gx_barcode_umi_groups.setdefault(genes[0], {}).setdefault(barcode, {}).setdefault(umi, []).append(read)
                if read.is_paired:
                    gx_barcode_umi_groups[genes[0]][barcode][umi].append(mate)

                if i % 10000 == 0:
                    # Call consensus for gene's whose bodies we've fully passed.
                    leftmost_start = start if not paired else next(iter(paired.values())).reference_start
                    for gene in list(gx_barcode_umi_groups.keys()):
                        gene_info = gene_infos[gene]
                        gene_contig = gene_info['chromosome']
                        gene_segment = gene_info['segment']
                        if (gene_contig < contig) or (gene_contig == contig and gene_segment.end <= leftmost_start):
                            barcode_umi_groups = gx_barcode_umi_groups.pop(gene)
                            for barcode, umi_groups in barcode_umi_groups.items():
                                for umi, reads in umi_groups.items():
                                    if len(reads) == 1:
                                        out.write(reads[0])
                                        continue

                                    tags = []
                                    if barcode_tag:
                                        tags.append((barcode_tag, barcode))
                                    if umi_tag:
                                        tags.append((umi_tag, umi))
                                    tags.extend([
                                        ('AS', sum(read.get_tag('AS') for read in reads)),
                                        ('NH', 1),
                                        ('HI', 1),
                                        ('GX', gene),
                                    ])
                                    gn = gene_info.get('gene_name')
                                    if gn:
                                        tags.append(('GN', gn))

                                    # Figure out what strand the consensus should map to
                                    consensus_strand = None
                                    gene_strand = gene_info['strand']
                                    if strand == 'forward':
                                        consensus_strand = gene_strand
                                    elif strand == 'reverse':
                                        consensus_strand = '-' if gene_strand == '+' else '+'

                                    # Save for multiprocessing later.
                                    args_q.put(([read.to_string()
                                                 for read in reads], header_dict, tags, consensus_strand))
                        else:
                            break

                    to_remove = 0
                    for gene in contig_gene_order[contig]:
                        if gene_infos[gene]['segment'].end <= leftmost_start:
                            to_remove += 1
                        else:
                            break
                    if to_remove > 0:
                        contig_gene_order[contig] = contig_gene_order[contig][to_remove:]

                    while True:
                        try:
                            result = results_q.get_nowait()
                            if result:
                                consensus = pysam.AlignedSegment.fromstring(result, header)
                                out.write(consensus)
                        except queue.Empty:
                            break

            for _ in range(len(workers)):
                args_q.put(None)

            for worker in workers:
                worker.join()

            while not results_q.empty():
                result = results_q.get()
                consensus = pysam.AlignedSegment.fromstring(result, header)
                out.write(consensus)
    # Sort and index
    return bam.sort_and_index_bam(temp_out_path, out_path, n_threads=n_threads, temp_dir=temp_dir)
