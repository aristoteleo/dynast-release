import array
import multiprocessing
import queue
from hashlib import sha256
from typing import Any, Dict, List, Optional

import ngs_tools as ngs
import numpy as np
import pysam
from tqdm import tqdm
from typing_extensions import Literal

from .. import config, utils
from ..logging import logger
from . import bam

BASES = ('A', 'C', 'G', 'T')
BASE_IDX = {base: i for i, base in enumerate(BASES)}


def call_consensus_from_reads(
    reads: List[pysam.AlignedSegment],
    header: pysam.AlignmentHeader,
    quality: int = 27,
    tags: Optional[Dict[str, Any]] = None,
    read_number: Optional[str] = None,
    shared_qname: Optional[str] = None
) -> pysam.AlignedSegment:
    """
    Call a single consensus alignment given a list of aligned reads.

    Reads must map to the same contig. Results are undefined otherwise.
    Additionally, consensus bases are called only for positions that match
    to the reference (i.e. no insertions allowed).

    If read_number is 'R1' or 'R2', we set the resulting consensus read
    to is_paired = True and is_read1/is_read2 = True accordingly,
    preserving the same QNAME for both consensus reads (if shared_qname is given).
    """
    if len(set(r.reference_name for r in reads)) > 1:
        raise Exception("Cannot call consensus from reads mapping to multiple contigs.")

    left_pos = min(r.reference_start for r in reads)
    right_pos = max(r.reference_end for r in reads)
    length = right_pos - left_pos

    sequence = np.zeros((length, len(BASES)), dtype=np.uint32)
    reference = np.full(length, -1, dtype=np.int8)  # -1 means unobserved

    for read in reads:
        read_seq = read.query_sequence.upper()
        read_quals = read.query_qualities
        for read_i, ref_i, ref_base in read.get_aligned_pairs(matches_only=False, with_seq=True):
            if ref_i is None or ref_base is None:
                continue
            i = ref_i - left_pos
            ref_base = ref_base.upper()
            if ref_base == 'N':
                continue
            if read_i is None:
                # Deletion
                if reference[i] < 0:
                    reference[i] = BASE_IDX[ref_base]
                continue
            read_base = read_seq[read_i]
            if read_base == 'N':
                continue
            if reference[i] < 0:
                reference[i] = BASE_IDX[ref_base]
            sequence[i, BASE_IDX[read_base]] += read_quals[read_i]

    # Build consensus
    consensus_mask = (sequence > 0).any(axis=1)
    consensus_length = consensus_mask.sum()
    consensus = np.zeros(consensus_length, dtype=np.uint8)
    qualities = np.zeros(consensus_length, dtype=np.uint8)

    cigar = []
    last_cigar_op = None
    cigar_n = 0
    md = []
    md_n = 0
    md_zero = True
    md_del = False
    nm = 0
    consensus_i = 0

    for i in range(length):
        ref_idx = reference[i]
        cigar_op = 'N'
        if ref_idx >= 0:
            seq = sequence[i]
            if (seq == 0).all():
                # Deletion
                cigar_op = 'D'
                if md_n > 0 or md_zero:
                    md.append(str(md_n))
                    md_n = 0
                if not md_del:
                    md.append('^')
                md.append(BASES[ref_idx])
                md_del = True
            else:
                md_del = False
                base_q = seq.max()
                if base_q < quality:
                    base = ref_idx
                else:
                    candidates = (seq == base_q).nonzero()[0]
                    if ref_idx in candidates:
                        base = ref_idx
                    else:
                        base = candidates[0]
                cigar_op = 'M'
                if base == ref_idx:
                    md_n += 1
                    md_zero = False
                else:
                    if md_n > 0 or md_zero:
                        md.append(str(md_n))
                        md_n = 0
                    md.append(BASES[ref_idx])
                    md_zero = True
                    nm += 1
                if consensus_mask[i]:
                    consensus[consensus_i] = base
                    qualities[consensus_i] = min(base_q, 42)
                    consensus_i += 1

        if cigar_op == last_cigar_op:
            cigar_n += 1
        else:
            if last_cigar_op:
                cigar.append(f"{cigar_n}{last_cigar_op}")
            last_cigar_op = cigar_op
            cigar_n = 1

    md.append(str(md_n))
    cigar.append(f"{cigar_n}{last_cigar_op}")

    al = pysam.AlignedSegment(header)
    if shared_qname is not None:
        al.query_name = shared_qname
    else:
        # default: hash together original QNAMEs
        all_names = ''.join(r.query_name for r in reads)
        al.query_name = sha256(all_names.encode('utf-8')).hexdigest()

    al.query_sequence = ''.join(BASES[b] for b in consensus)
    al.query_qualities = array.array('B', qualities)
    al.reference_name = reads[0].reference_name
    al.reference_id = reads[0].reference_id
    al.reference_start = left_pos
    al.mapping_quality = 255
    al.cigarstring = ''.join(cigar)

    # Add tags
    tags = tags or {}
    tags.update({'MD': ''.join(md), 'NM': nm})
    al.set_tags(list(tags.items()))

    # Mark R1 or R2 if provided
    if read_number == 'R1':
        al.is_paired = True
        al.is_read1 = True
        al.is_read2 = False
    elif read_number == 'R2':
        al.is_paired = True
        al.is_read1 = False
        al.is_read2 = True
    else:
        al.is_paired = False

    # Force false for these
    al.is_unmapped = False
    al.is_duplicate = False
    al.is_qcfail = False
    al.is_secondary = False
    al.is_supplementary = False

    return al


def call_consensus_from_reads_process(
    reads,
    header,
    tags,
    strand=None,
    read_number=None,
    shared_qname=None,
    quality=27
):
    """
    Helper for multiprocessing calls.
    """
    header = pysam.AlignmentHeader.from_dict(header)
    pysam_reads = [pysam.AlignedSegment.fromstring(r, header) for r in reads]

    aln = call_consensus_from_reads(
        pysam_reads,
        header,
        quality=quality,
        tags=tags,
        read_number=read_number,
        shared_qname=shared_qname
    )
    if strand == '-':
        aln.is_reverse = True
    return aln.to_string()


def consensus_worker(args_q, results_q, quality=27):
    """
    Worker that reads tasks from args_q, calls call_consensus_from_reads_process.
    """
    while True:
        try:
            _args = args_q.get(timeout=1)
        except queue.Empty:
            continue
        if _args is None:
            return
        results_q.put(call_consensus_from_reads_process(*_args, quality=quality))


def get_read_number(read: pysam.AlignedSegment) -> str:
    """
    Return 'R1' if read.is_read1 is True, else 'R2'.
    For single-end reads, default to 'R1'.
    """
    if read.is_paired:
        return 'R1' if read.is_read1 else 'R2'
    else:
        return 'R1'


def call_consensus(
    bam_path: str,
    out_path: str,
    gene_infos: dict,
    strand: Literal['forward', 'reverse', 'unstranded'] = 'forward',
    umi_tag: Optional[str] = None,
    barcode_tag: Optional[str] = None,
    gene_tag: str = 'GX',
    barcodes: Optional[List[str]] = None,
    quality: int = 27,
    add_RS_RI: bool = False,
    temp_dir: Optional[str] = None,
    n_threads: int = 8,
    collapse_r1_r2: bool = False  # <-- ### ADDED
) -> str:
    """
    Call consensus sequences from BAM.

    If collapse_r1_r2 is True, then R1 and R2 from the same UMI
    are stored together in one group, producing a single consensus.

    If collapse_r1_r2 is False, R1 and R2 from the same UMI
    are stored separately, producing two consensus reads (one for R1, one for R2).

    Args:
        bam_path: Path to BAM
        out_path: Output BAM path
        gene_infos: Gene info from GTF
        strand: Protocol strandedness
        umi_tag: BAM tag for UMI
        barcode_tag: BAM tag for barcode
        gene_tag: BAM tag for assigned gene
        barcodes: optional filter for barcodes
        quality: consensus base quality threshold
        add_RS_RI: optional debug tags
        temp_dir: optional temp dir
        n_threads: number of threads
        collapse_r1_r2: if True, combine R1/R2. if False, produce separate consensus.
    """

    def skip_alignment(read, required):
        return (
            read.is_secondary
            or read.is_unmapped
            or any(not read.has_tag(t) for t in required)
        )

    def find_genes(contig, start, end, read_strand=None):
        genes = []
        for gene in contig_gene_order.get(contig, []):
            if read_strand and read_strand != gene_infos[gene]['strand']:
                continue
            gene_segment = gene_infos[gene]['segment']
            if end <= gene_segment.start:
                break
            if start >= gene_segment.start and end <= gene_segment.end:
                genes.append(gene)
        return genes

    def swap_gene_tags(r, gene):
        t = dict(r.get_tags())
        if gene_tag and r.has_tag(gene_tag):
            del t[gene_tag]
        t['GX'] = gene
        gn = gene_infos.get(gene, {}).get('gene_name')
        if gn:
            t['GN'] = gn
        r.set_tags(list(t.items()))
        return r

    def create_tags_and_strand(barcode, umi, reads, ginfo):
        as_sum = sum(r.get_tag('AS') for r in reads if r.has_tag('AS'))
        tags = {
            'AS': as_sum,
            'NH': 1,
            'HI': 1,
            config.BAM_CONSENSUS_READ_COUNT_TAG: len(reads),
        }
        if barcode_tag:
            tags[barcode_tag] = barcode
        if umi_tag:
            tags[umi_tag] = umi

        if gene_tag:
            gene_id = None
            for rr in reads:
                if rr.has_tag(gene_tag):
                    gene_id = rr.get_tag(gene_tag)
                    break
            if gene_id:
                tags['GX'] = gene_id
                gn = ginfo.get('gene_name')
                if gn:
                    tags['GN'] = gn

        if add_RS_RI:
            tags['RS'] = ';'.join(r.query_name for r in reads)
            tags['RI'] = ';'.join(
                str(r.get_tag('HI')) if r.has_tag('HI') else '0' for r in reads
            )

        consensus_strand = None
        gene_strand = ginfo['strand']
        if strand == 'forward':
            consensus_strand = gene_strand
        elif strand == 'reverse':
            consensus_strand = '-' if gene_strand == '+' else '+'

        return tags, consensus_strand

    if add_RS_RI:
        logger.warning("RS and RI tags may greatly increase BAM size.")

    # Build index of genes per contig
    contig_gene_order = {}
    for gene_id, ginfo in gene_infos.items():
        contig_gene_order.setdefault(ginfo['chromosome'], []).append(gene_id)
    for contig in list(contig_gene_order.keys()):
        contig_gene_order[contig] = sorted(
            contig_gene_order[contig],
            key=lambda g: tuple(gene_infos[g]['segment'])
        )

    # ### CHANGED
    # Instead of storing R1 and R2 in separate keys unconditionally,
    # we let the user decide via 'collapse_r1_r2'.
    #
    # Data structure:
    #   gx_barcode_umi_groups[gene][barcode][umi][subkey] -> list_of_reads
    # where 'subkey' is either:
    #   'ALL'   (if collapse_r1_r2=True)
    #   'R1' or 'R2' (if collapse_r1_r2=False)
    gx_barcode_umi_groups = {}
    paired = {}

    required_tags = []
    if umi_tag:
        required_tags.append(umi_tag)
    if barcode_tag:
        required_tags.append(barcode_tag)

    manager = multiprocessing.Manager()
    args_q = manager.Queue(1000 * n_threads)
    results_q = manager.Queue()

    workers = [
        multiprocessing.Process(
            target=consensus_worker,
            args=(args_q, results_q, quality),
            daemon=True
        )
        for _ in range(n_threads)
    ]
    for w in workers:
        w.start()

    temp_out_path = utils.mkstemp(dir=temp_dir)
    with pysam.AlignmentFile(bam_path, 'rb') as f_in:
        header_dict = f_in.header.to_dict()
        hd = header_dict.setdefault('HD', {'VN': '1.4', 'SO': 'unsorted'})
        hd['SO'] = 'unsorted'
        header = pysam.AlignmentHeader.from_dict(header_dict)

        total_reads = ngs.bam.count_bam(bam_path)
        with pysam.AlignmentFile(temp_out_path, 'wb', header=header) as out:
            for i, read in tqdm(enumerate(f_in.fetch()),
                                total=total_reads, ascii=True, smoothing=0.01,
                                desc="Calling consensus"):
                if skip_alignment(read, required_tags):
                    continue

                barcode = read.get_tag(barcode_tag) if barcode_tag else None
                if barcode == '-' or (barcodes and barcode not in barcodes):
                    continue

                umi = read.get_tag(umi_tag) if umi_tag else None
                read_id = read.query_name
                alignment_index = read.get_tag('HI') if read.has_tag('HI') else 1
                key = (read_id, alignment_index)

                # Handle pairing
                mate = None
                if read.is_paired:
                    if key not in paired:
                        paired[key] = read
                        continue
                    mate = paired.pop(key)
                    if not umi:
                        umi = (mate.reference_start, mate.reference_end)

                # Determine the gene
                start = read.reference_start
                end = read.reference_end
                read_strand = None
                if strand in ('forward', 'reverse'):
                    if read.is_paired:
                        if strand == 'forward':
                            read_strand = '+' if not read.is_reverse else '-'
                        else:  # 'reverse'
                            read_strand = '-' if not read.is_reverse else '+'
                    else:
                        if strand == 'forward':
                            read_strand = '-' if read.is_reverse else '+'
                        else:
                            read_strand = '+' if read.is_reverse else '-'

                if gene_tag and read.has_tag(gene_tag):
                    genes = [read.get_tag(gene_tag)]
                else:
                    contig = read.reference_name
                    genes = find_genes(contig, start, end, read_strand)

                if len(genes) != 1:
                    # If not exactly one gene, write raw
                    out.write(read)
                    if mate:
                        out.write(mate)
                    continue

                gene = genes[0]

                # ### CHANGED
                # Decide subkey = 'ALL' (if collapsing) or 'R1'/'R2' (if separating)
                if collapse_r1_r2:
                    subkey = 'ALL'
                else:
                    subkey = get_read_number(read)

                gx_barcode_umi_groups \
                    .setdefault(gene, {}) \
                    .setdefault(barcode, {}) \
                    .setdefault(umi, {}) \
                    .setdefault(subkey, []) \
                    .append(read)

                # Also store mate if present
                if mate:
                    if collapse_r1_r2:
                        mate_subkey = 'ALL'
                    else:
                        mate_subkey = get_read_number(mate)
                    gx_barcode_umi_groups[gene][barcode][umi] \
                        .setdefault(mate_subkey, []) \
                        .append(mate)

                # Periodically flush old genes
                if i % 10000 == 0:
                    leftmost = min(
                        read.reference_start,
                        mate.reference_start if mate else read.reference_start
                    )
                    for g in list(gx_barcode_umi_groups.keys()):
                        ginfo = gene_infos[g]
                        gene_contig = ginfo['chromosome']
                        gene_segment = ginfo['segment']
                        if (gene_contig < read.reference_name) or (
                            gene_contig == read.reference_name and gene_segment.end <= leftmost
                        ):
                            bc_map = gx_barcode_umi_groups.pop(g)
                            for bc, umi_map in bc_map.items():
                                for this_umi, sub_map in umi_map.items():
                                    # Build a stable QNAME for everything in sub_map
                                    all_names = []
                                    for subkey_, reads_list in sub_map.items():
                                        all_names.extend(r.query_name for r in reads_list)
                                    shared_qname = sha256(''.join(all_names).encode('utf-8')).hexdigest()

                                    for subk, reads_list in sub_map.items():
                                        if len(reads_list) == 1:
                                            single_read = swap_gene_tags(reads_list[0], g)
                                            single_read.query_name = shared_qname
                                            if not collapse_r1_r2:
                                                # If not collapsing, check if subk = 'R1'/'R2'
                                                if subk == 'R1':
                                                    single_read.is_paired = True
                                                    single_read.is_read1 = True
                                                    single_read.is_read2 = False
                                                elif subk == 'R2':
                                                    single_read.is_paired = True
                                                    single_read.is_read1 = False
                                                    single_read.is_read2 = True
                                            out.write(single_read)
                                        else:
                                            tags, cstrand = create_tags_and_strand(bc, this_umi, reads_list, ginfo)
                                            # subk might be 'ALL' or 'R1'/'R2'
                                            args_q.put((
                                                [r.to_string() for r in reads_list],
                                                header_dict,
                                                tags,
                                                cstrand,
                                                subk if not collapse_r1_r2 else None,  # read_number if separate
                                                shared_qname
                                            ))

                    # Drain queue
                    while True:
                        try:
                            result = results_q.get_nowait()
                            if result:
                                out.write(pysam.AlignedSegment.fromstring(result, header))
                        except queue.Empty:
                            break

            # Final flush
            for g, bc_map in gx_barcode_umi_groups.items():
                ginfo = gene_infos[g]
                for bc, umi_map in bc_map.items():
                    for this_umi, sub_map in umi_map.items():
                        all_names = []
                        for subk, reads_list in sub_map.items():
                            all_names.extend(r.query_name for r in reads_list)
                        shared_qname = sha256(''.join(all_names).encode('utf-8')).hexdigest()

                        for subk, reads_list in sub_map.items():
                            if len(reads_list) == 1:
                                single_read = swap_gene_tags(reads_list[0], g)
                                single_read.query_name = shared_qname
                                if not collapse_r1_r2:
                                    if subk == 'R1':
                                        single_read.is_paired = True
                                        single_read.is_read1 = True
                                        single_read.is_read2 = False
                                    elif subk == 'R2':
                                        single_read.is_paired = True
                                        single_read.is_read1 = False
                                        single_read.is_read2 = True
                                out.write(single_read)
                            else:
                                tags, cstrand = create_tags_and_strand(bc, this_umi, reads_list, ginfo)
                                args_q.put((
                                    [r.to_string() for r in reads_list],
                                    header_dict,
                                    tags,
                                    cstrand,
                                    subk if not collapse_r1_r2 else None,
                                    shared_qname
                                ))

            # Signal termination
            for _ in range(len(workers)):
                args_q.put(None)
            for w in workers:
                w.join()

            # Gather last results
            while not results_q.empty():
                result = results_q.get()
                out.write(pysam.AlignedSegment.fromstring(result, header))

    # Sort and index as usual
    return bam.sort_and_index_bam(temp_out_path, out_path, n_threads=n_threads, temp_dir=temp_dir)
