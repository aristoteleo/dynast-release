import argparse
import logging
import os
import shutil
import sys
import warnings
from typing import Optional

from . import __version__
from .config import BAM_GENE_TAG
from .logging import logger
from .preprocessing import CONVERSION_COMPLEMENT
from .technology import STRAND_MAP, TECHNOLOGIES_MAP, detect_strand
from .utils import flatten_iter, patch_mp_connection_bpo_17560


def print_technologies():
    """Displays a list of supported technologies along with whether a whitelist
    is provided for that technology.
    """
    headers = ['name', 'whitelist', 'barcode', 'umi', 'cDNA', 'strand']
    rows = [headers]

    print('List of supported single-cell technologies\n')
    print('Positions syntax: `input file index, start position, end position`')
    print('When start & end positions are None, refers to the entire file\n')
    for key in sorted(TECHNOLOGIES_MAP):
        t = TECHNOLOGIES_MAP[key]
        chem = t.chemistry
        row = [
            t.name, 'yes' if chem.has_whitelist else '',
            ' '.join(str(_def) for _def in chem.cell_barcode_parser) if chem.has_cell_barcode else '',
            ' '.join(str(_def) for _def in chem.umi_parser) if chem.has_umi else '',
            ' '.join(str(_def) for _def in chem.cdna_parser), chem.strand.name
        ]
        rows.append(row)

    max_lens = []
    for i in range(len(headers)):
        max_lens.append(len(headers[i]))
        for row in rows[1:]:
            max_lens[i] = max(max_lens[i], len(row[i]))

    rows.insert(1, ['-' * l for l in max_lens])  # noqa
    for row in rows:
        for col, l in zip(row, max_lens):
            print(col.ljust(l + 4), end='')
        print()
    sys.exit(1)


def setup_ref_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `ref` command.

    Args:
        parser: Argparse parser to add the `ref` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_ref = parser.add_parser(
        'ref',
        description='Build a STAR index from a reference',
        help='Build a STAR index from a reference',
        parents=[parent],
    )
    parser_ref._actions[0].help = parser_ref._actions[0].help.capitalize()

    required_ref = parser_ref.add_argument_group('required arguments')
    required_ref.add_argument(
        '-i',
        metavar='INDEX',
        help='Path to the directory where the STAR index will be generated',
        type=str,
        required=True
    )
    parser_ref.add_argument(
        '-m', metavar='MEMORY', help='Maximum memory used, in GB (default: 16)', type=int, default=16
    )
    parser_ref.add_argument(
        'fasta',
        help='Genomic FASTA file',
        type=str,
    )
    parser_ref.add_argument(
        'gtf',
        help='Reference GTF file',
        type=str,
    )

    return parser_ref


def setup_align_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `align` command.

    Args:
        parser: Argparse parser to add the `align` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_align = parser.add_parser(
        'align',
        description='Align FASTQs',
        help='Align FASTQs',
        parents=[parent],
    )
    parser_align._actions[0].help = parser_align._actions[0].help.capitalize()

    required_align = parser_align.add_argument_group('required arguments')
    required_align.add_argument(
        '-i', metavar='INDEX', help='Path to the directory where the STAR index is located', type=str, required=True
    )
    required_align.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    required_align.add_argument(
        '-x',
        metavar='TECHNOLOGY',
        help='Single-cell technology used. `dynast --list` to view all supported technologies',
        type=str,
        required=True,
        choices=TECHNOLOGIES_MAP.keys()
    )
    parser_align.add_argument(
        '--strand',
        help='Read strandedness. By default, this is auto-detected depending on the technology.',
        choices=['forward', 'reverse', 'unstranded', 'auto'],
        default='auto',
    )
    parser_align.add_argument(
        '-w',
        metavar='WHITELIST',
        help=('Path to file of whitelisted barcodes to correct to. '
              'If not provided, all barcodes are used.'),
    )
    parser_align.add_argument('--overwrite', help='Overwrite existing files.', action='store_true')
    parser_align.add_argument(
        '--STAR-overrides', metavar='ARGUMENTS', help='Arguments to pass directly to STAR.', type=str, default=None
    )
    parser_align.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )
    parser_align.add_argument(
        'fastqs',
        help=(
            'FASTQ files, where the first contains biological reads and the second contains barcode and UMI reads. '
            'If `-x smartseq`, this is a single manifest CSV file where '
            'the first column contains cell IDs and the next two columns contain '
            'paths to FASTQs (the third column may contain a dash `-` for single-end reads).'
        ),
        nargs='+'
    )

    return parser_align


def setup_consensus_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `consensus` command.

    Args:
        parser: Argparse parser to add the `consensus` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_consensus = parser.add_parser(
        'consensus',
        description='Generate consensus sequences',
        help='Generate consensus sequences',
        parents=[parent],
    )
    parser_consensus._actions[0].help = parser_consensus._actions[0].help.capitalize()

    required_consensus = parser_consensus.add_argument_group('required arguments')
    required_consensus.add_argument(
        '-g', metavar='GTF', help='Path to GTF file used to generate the STAR index', type=str, required=True
    )
    parser_consensus.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_consensus.add_argument(
        '--umi-tag',
        metavar='TAG',
        help=(
            'BAM tag to use as unique molecular identifiers (UMI). If not provided, '
            'all reads are assumed to be unique. (default: None)'
        ),
        type=str,
        default=None,
    )
    parser_consensus.add_argument(
        '--barcode-tag',
        metavar='TAG',
        help=(
            'BAM tag to use as cell barcodes. If not provided, all reads are '
            'assumed to be from a single cell. (default: None)'
        ),
        type=str,
        default=None,
    )
    parser_consensus.add_argument(
        '--gene-tag',
        metavar='TAG',
        help=f'BAM tag to use as gene assignments (default: {BAM_GENE_TAG})',
        type=str,
        default=BAM_GENE_TAG,
    )
    parser_consensus.add_argument(
        '--strand',
        help='Read strandedness. By default, this is auto-detected from the BAM.',
        choices=['forward', 'reverse', 'unstranded', 'auto'],
        default='auto',
    )
    parser_consensus.add_argument(
        '--quality',
        metavar='QUALITY',
        help=(
            'Base quality threshold. When generating a consensus nucleotide at a '
            'certain position, the base with smallest error probability below this '
            'quality threshold is chosen. If no base meets this criteria, the '
            'reference base is chosen. (default: 27)'
        ),
        type=int,
        default=27
    )
    parser_consensus.add_argument(
        '--barcodes',
        metavar='TXT',
        help=('Textfile containing filtered cell barcodes. Only these barcodes will '
              'be processed.'),
        type=str,
    )
    parser_consensus.add_argument(
        '--add-RS-RI',
        help=(
            'Add custom RS and RI tags to the output BAM, each of which contain a '
            'semi-colon delimited list of read names (RS) and alignment indices (RI) '
            'of the reads and alignments from which the consensus is derived. This '
            'option is useful for debugging.'
        ),
        action='store_true'
    )
    parser_consensus.add_argument(
        'bam',
        help=(
            'Alignment BAM file that contains the appropriate UMI and barcode tags, '
            'specifiable with `--umi-tag`, and `--barcode-tag`.'
        ),
        type=str,
    )
    return parser_consensus


def setup_count_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `count` command.

    Args:
        parser: Argparse parser to add the `count` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_count = parser.add_parser(
        'count',
        description='Quantify unlabeled and labeled RNA',
        help='Quantify unlabeled and labeled RNA',
        parents=[parent],
    )
    parser_count._actions[0].help = parser_count._actions[0].help.capitalize()

    required_count = parser_count.add_argument_group('required arguments')
    required_count.add_argument(
        '-g', metavar='GTF', help='Path to GTF file used to generate the STAR index', type=str, required=True
    )
    required_count.add_argument(
        '--conversion',
        metavar='CONVERSION',
        help=(
            'The type of conversion(s) introduced at a single timepoint. Multiple conversions '
            'can be specified with a comma-delimited list. For example, T>C and A>G is TC,AG. '
            'This option can be specified multiple times (i.e. dual labeling), for each labeling '
            'timepoint.'
        ),
        action='append',
        required=True,
        choices=sorted(CONVERSION_COMPLEMENT.keys())
    )
    parser_count.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_count.add_argument(
        '--umi-tag',
        metavar='TAG',
        help=(
            'BAM tag to use as unique molecular identifiers (UMI). If not provided, '
            'all reads are assumed to be unique. (default: None)'
        ),
        type=str,
        default=None,
    )
    parser_count.add_argument(
        '--barcode-tag',
        metavar='TAG',
        help=(
            'BAM tag to use as cell barcodes. If not provided, all reads are '
            'assumed to be from a single cell. (default: None)'
        ),
        type=str,
        default=None,
    )
    parser_count.add_argument(
        '--gene-tag',
        metavar='TAG',
        help=f'BAM tag to use as gene assignments (default: {BAM_GENE_TAG})',
        type=str,
        default=BAM_GENE_TAG,
    )
    parser_count.add_argument(
        '--strand',
        help='Read strandedness. By default, this is auto-detected from the BAM.',
        choices=['forward', 'reverse', 'unstranded', 'auto'],
        default='auto',
    )
    parser_count.add_argument(
        '--quality',
        metavar='QUALITY',
        help=(
            'Base quality threshold. Only bases with PHRED quality greater than '
            'this value will be considered when counting conversions. (default: 27)'
        ),
        type=int,
        default=27
    )
    parser_count.add_argument(
        '--snp-threshold',
        metavar='THRESHOLD',
        help=(
            'Conversions with (# conversions) / (# reads) greater than this '
            'threshold will be considered a SNP and ignored. '
            '(default: no SNP detection)'
        ),
        type=float,
        default=None,
    )
    parser_count.add_argument(
        '--snp-min-coverage',
        metavar='THRESHOLD',
        help=(
            'For a conversion to be considered as a SNP, there must be at least '
            'this many reads mapping to that region. (default: 1)'
        ),
        type=int,
        default=1,
    )
    parser_count.add_argument(
        '--snp-csv',
        metavar='CSV',
        help=('CSV file of two columns: contig (i.e. chromosome) and genome position '
              'of known SNPs'),
        type=str,
        default=None,
    )
    parser_count.add_argument(
        '--barcodes',
        metavar='TXT',
        help=('Textfile containing filtered cell barcodes. Only these barcodes will '
              'be processed.'),
        type=str,
    )
    parser_count.add_argument(
        '--gene-names',
        help=('Group counts by gene names instead of gene IDs when generating the h5ad file.'),
        action='store_true'
    )

    splicing_group = parser_count.add_mutually_exclusive_group()
    splicing_group.add_argument(
        '--no-splicing',
        '--transcriptome-only',
        help=(
            'Do not assign reads a splicing status (spliced, unspliced, ambiguous) '
            'and ignore reads that are not assigned to the transcriptome.'
        ),
        action='store_true'
    )
    splicing_group.add_argument(
        '--exon-overlap',
        help=(
            'Algorithm to use to detect spliced reads (that overlap exons). '
            'May be `strict`, which assigns reads as spliced if it only overlaps '
            'exons, or `lenient`, which assigns reads as spliced if it does not '
            'overlap with any introns of at least one transcript. (default: strict)'
        ),
        type=str,
        choices=['lenient', 'strict'],
        default='strict',
    )
    parser_count.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )
    parser_count.add_argument(
        '--control',
        help='Indicate this is a control sample, which is used to detect SNPs.',
        action='store_true',
    )
    parser_count.add_argument(
        '--dedup-mode',
        help=(
            'Deduplication mode for UMI-based technologies (required `--umi-tag`). '
            'Available choices are: `auto`, `conversion`, `exon`. When `conversion` is used, '
            'reads that have at least one of the provided conversions is prioritized. '
            'When `exon` is used, exonic reads are prioritized. By default (`auto`), '
            'the BAM is inspected to select the appropriate mode.'
        ),
        type=str,
        choices=['auto', 'conversion', 'exon'],
    )
    parser_count.add_argument(
        '--overwrite',
        help='Overwrite existing files.',
        action='store_true',
    )
    parser_count.add_argument(
        'bam',
        help=(
            'Alignment BAM file that contains the appropriate UMI and barcode tags, '
            'specifiable with `--umi-tag`, and `--barcode-tag`.'
        ),
        type=str,
    )

    return parser_count


def setup_estimate_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `estimate` command.

    Args:
        parser: Argparse parser to add the `estimate` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_estimate = parser.add_parser(
        'estimate',
        description='Estimate fraction of labeled RNA',
        help='Estimate fraction of labeled RNA',
        parents=[parent],
    )
    parser_estimate._actions[0].help = parser_estimate._actions[0].help.capitalize()

    parser_estimate.add_argument(
        '--reads',
        help=(
            'Read groups to perform estimation on. '
            'This option can be used multiple times to estimate multiple groups. '
            '(default: all possible reads groups)'
        ),
        action='append',
        choices=['total', 'transcriptome', 'spliced', 'unspliced'],
        default=None
    )

    parser_estimate.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_estimate.add_argument(
        '--groups',
        metavar='CSV',
        help=(
            'CSV containing cell (barcode) groups, where the first column is the barcode '
            'and the second is the group name the cell belongs to. Cells will be combined per '
            'group for estimation of parameters specified by `--groups-for`.'
        ),
        action='append',
        default=None,
    )
    parser_estimate.add_argument(
        '--ignore-groups-for-pi',
        help=(
            'Ignore cell groupings when estimating the fraction of labeled RNA. '
            'This option only has an effect when `--groups` is also specified.'
        ),
        action='store_true',
    )
    parser_estimate.add_argument(
        '--genes',
        metavar='TXT',
        help=('Textfile containing list of genes to use. All other genes will be '
              'treated as if they do not exist.'),
        type=str,
    )
    parser_estimate.add_argument(
        '--cell-threshold',
        metavar='COUNT',
        help='A cell must have at least this many reads for correction. (default: 1000)',
        type=int,
        default=1000,
    )
    parser_estimate.add_argument(
        '--cell-gene-threshold',
        metavar='COUNT',
        help='A cell-gene pair must have at least this many reads for correction. (default: 16)',
        type=int,
        default=16
    )
    parser_estimate.add_argument(
        '--gene-names',
        help=('Group counts by gene names instead of gene IDs when generating h5ad file'),
        action='store_true'
    )
    parser_estimate.add_argument(
        '--downsample',
        metavar='NUM',
        help=(
            'Downsample the number of reads (UMIs). If a decimal between 0 and 1 is given, '
            'then the number is interpreted as the proportion of remaining reads. If an integer '
            'is given, the number is interpreted as the absolute number of remaining reads.'
        ),
        type=float,
        default=None
    )
    parser_estimate.add_argument(
        '--downsample-mode',
        metavar='MODE',
        help=(
            'Downsampling mode. Can be one of: `uniform`, `cell`, `group`. If `uniform`, all reads '
            '(UMIs) are downsampled uniformly at random. If `cell`, only cells that have more '
            'reads than the argument to `--downsample` are downsampled to exactly that number. '
            'If `group`, identical to `cell` but per group specified by `--groups`.'
        ),
        type=str,
        choices=['uniform', 'cell', 'group'],
        default='uniform'
    )
    parser_estimate.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )
    parser_estimate.add_argument('--seed', help=argparse.SUPPRESS, type=int, default=None)
    parser_estimate.add_argument(
        '--control',
        help=('Indicate this is a control sample, only the background mutation rate '
              'will be estimated.'),
        action='store_true',
    )
    parser_estimate.add_argument(
        '--p-e',
        help='Textfile containing a single number, indicating the estimated background mutation rate',
        type=str,
        default=None,
    )
    parser_estimate.add_argument(
        'count_dirs',
        help=(
            'Path to directory that contains `dynast count` output. When multiple are provided, '
            'the barcodes in each of the count directories are suffixed with `-i` where i is '
            'a 0-indexed integer.'
        ),
        type=str,
        nargs='+',
    )

    return parser_estimate


def parse_ref(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `ref` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    if os.path.exists(args.i):
        parser.error(
            f'STAR index directory {args.i} already exists. '
            'Please provide a different directory or remove the existing one.'
        )

    from .ref import ref
    ref(
        args.fasta,
        args.gtf,
        args.i,
        n_threads=args.t,
        memory=args.m * (1024**3),
        temp_dir=temp_dir,
    )


def parse_align(parser, args, temp_dir=None):
    technology = TECHNOLOGIES_MAP[args.x]
    if technology.name == 'smartseq':
        if len(args.fastqs) != 1:
            parser.error('A single manifest TSV must be provided for technology `smartseq`')
        with open(args.fastqs[0], 'r') as f:
            cell_ids = []
            for line in f:
                if line.isspace():
                    continue

                cell_id, fastq_1, fastq_2 = line.strip().split(',')
                if cell_id in cell_ids:
                    parser.error(f'Found duplicate cell ID in manifest CSV: {cell_id}')
                if not os.path.exists(fastq_1):
                    parser.error(
                        f'{fastq_1} does not exist. '
                        f'All paths in {args.fastqs[0]} must be relative to the '
                        'current working directory or absolute.'
                    )
                if fastq_2 != '-' and not os.path.exists(fastq_2):
                    parser.error(
                        f'{fastq_2} does not exist. '
                        f'All paths in {args.fastqs[0]} must be relative to the '
                        'current working directory or absolute.'
                    )
                cell_ids.append(cell_id)
    else:
        # Check number of fastqs
        if len(args.fastqs) != 2:
            parser.error(f'Two input FASTQs were expected, but {len(args.fastqs)} were provided')

    # STAR overrides
    overrides = {}
    if args.STAR_overrides:
        # Strip any quotes or double-quotes because those shouldn't be part of the string
        args.STAR_overrides = args.STAR_overrides.strip('\"\'')
        arg = None
        for part in args.STAR_overrides.split(' '):
            if part.startswith(('-', '--')):
                # In case it is a number
                try:
                    float(part)
                    overrides.setdefault(arg, []).append(part)
                    continue
                except ValueError:
                    pass
                arg = part
            else:
                overrides.setdefault(arg, []).append(part)
        # Clean
        overrides = {arg: parts[0] if len(parts) == 1 else parts for arg, parts in overrides.items()}

    # Detect strand
    strand = args.strand
    if args.strand == 'auto':
        strand = STRAND_MAP[technology.chemistry.strand]
        logger.info(
            f'Auto-detected strandedness `{strand}` for technology `{technology.name}`. '
            'Use `--strand` to override.'
        )

    from .align import align
    align(
        args.fastqs,
        args.i,
        args.o,
        technology,
        whitelist_path=args.w,
        strand=strand,
        n_threads=args.t,
        temp_dir=temp_dir,
        nasc=args.nasc,
        overrides=overrides
    )


def parse_consensus(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `consensus` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    # Check quality
    if args.quality < 0 or args.quality > 41:
        parser.error('`--quality` must be in [0, 42)')

    # Read barcodes
    barcodes = set()
    if args.barcodes:
        if not args.barcode_tag:
            parser.error('`--barcodes` may only be provided with `--barcode-tag`.')

        with open(args.barcodes, 'r') as f:
            barcodes = set(line.strip() for line in f if not line.isspace())
        logger.warning(f'Ignoring cell barcodes not in the {len(barcodes)} barcodes provided by `--barcodes`')
    else:
        logger.warning('`--barcodes` not provided. All cell barcodes will be processed.')

    # Detect strand
    strand = args.strand
    if args.strand == 'auto':
        strand = detect_strand(args.bam)
        if strand is None:
            parser.error('Failed to auto-detect strandedness from BAM. Use `--strand` to provide strandedness.')
        else:
            logger.info(f'Auto-detected strandedness: {strand}. Use `--strand` to override.')

    from .consensus import consensus
    consensus(
        args.bam,
        args.g,
        args.o,
        strand=strand,
        umi_tag=args.umi_tag,
        barcode_tag=args.barcode_tag,
        gene_tag=args.gene_tag,
        barcodes=barcodes,
        quality=args.quality,
        add_RS_RI=args.add_RS_RI,
        n_threads=args.t,
        temp_dir=temp_dir,
    )


def parse_count(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `count` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    # Check quality
    if args.quality < 0 or args.quality > 41:
        parser.error('`--quality` must be in [0, 42)')

    # Read barcodes
    barcodes = set()
    if args.barcodes:
        if not args.barcode_tag:
            parser.error('`--barcodes` may only be provided with `--barcode-tag`.')

        with open(args.barcodes, 'r') as f:
            barcodes = set(line.strip() for line in f if not line.isspace())
        logger.warning(f'Ignoring cell barcodes not in the {len(barcodes)} barcodes provided by `--barcodes`')
    else:
        logger.warning('`--barcodes` not provided. All cell barcodes will be processed.')

    if not args.gene_tag:
        parser.error('`--gene-tag` must be a valid tag')

    # Check conversions
    if args.nasc and args.conversion != ['TC']:
        parser.error('`--conversion TC` is required for `--nasc`')
    conversions = []
    for arg in args.conversion:
        conversions.append(arg.upper().split(','))
    flattened = list(flatten_iter(conversions))
    if len(set(flattened)) != len(flattened):
        parser.error('duplicate conversions are not allowed for `--conversion`')

    # Convert conversions to frozenset of tuples.
    conversions = frozenset(frozenset(conv) for conv in conversions)

    # Detect strand
    strand = args.strand
    if args.strand == 'auto':
        strand = detect_strand(args.bam)
        if strand is None:
            parser.error('Failed to auto-detect strandedness from BAM. Use `--strand` to provide strandedness.')
        else:
            logger.info(f'Auto-detected strandedness: {strand}. Use `--strand` to override.')

    from .count import count
    count(
        args.bam,
        args.g,
        args.o,
        strand=strand,
        umi_tag=args.umi_tag,
        barcode_tag=args.barcode_tag,
        gene_tag=args.gene_tag,
        barcodes=barcodes,
        control=args.control,
        quality=args.quality,
        conversions=conversions,
        snp_threshold=args.snp_threshold,
        snp_min_coverage=args.snp_min_coverage,
        snp_csv=args.snp_csv,
        n_threads=args.t,
        temp_dir=temp_dir,
        nasc=args.nasc,
        overwrite=args.overwrite,
        velocity=not args.no_splicing,
        strict_exon_overlap=args.exon_overlap == 'strict',
        dedup_mode=args.dedup_mode,
        by_name=args.gene_names
    )


def parse_estimate(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `estimate` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    # Check control constraints
    if args.control and args.p_e:
        parser.error('`--control` and `--p-e` can not be used together')

    # Check p_e is in correct format (only a single number)
    control_p_e = None
    if args.p_e:
        with open(args.p_e, 'r') as f:
            try:
                control_p_e = float(f.read().strip())
            except ValueError:
                parser.error('`--p-e` must be a textfile containing a single decimal number')

    # group CSV must be proivded when there are multiple input directories
    if len(args.count_dirs) > 1 and not args.groups:
        parser.error('`--group` CSVs must be provided when using multiple input directories')

    # If group CSV(s) are provided, the number must match input directories
    if args.groups and len(args.groups) != len(args.count_dirs):
        parser.error('Number of `--group` CSVs must match number of input directories')

    # Multiple count dirs can't be used with nasc
    if len(args.count_dirs) > 1 and args.nasc:
        parser.error('`--nasc` does not support multiple input directories')

    # Read genes
    genes = None
    if args.genes:
        if args.by_name:
            logger.warning(
                '`--genes` were provided with `--gene-names`. '
                'Make sure your gene list contains gene names instead of IDs. '
                'IDs should be used for any genes that do not have a name.'
            )
        with open(args.genes, 'r') as f:
            genes = [line.strip() for line in f if not line.isspace()]
        logger.warning(f'Ignoring genes not in the {len(genes)} genes provided by `--genes`')

    # War

    # Parse cell groups csv(s)
    groups = []
    if args.groups:
        for path in args.groups:
            groups_part = {}
            with open(path, 'r') as f:
                for line in f:
                    if line.isspace():
                        continue
                    barcode, group = line.strip().split(',')

                    if barcode in groups_part:
                        parser.error(f'Found duplicate barcode {barcode} in {path}')

                    groups_part[barcode] = group
            groups.append(groups_part)

    if args.ignore_groups_for_pi and not (args.groups):
        parser.error('`--ignore-groups-for-pi` can not be used without `--groups`')

    if not args.reads:
        args.reads = 'complete'

    # Check that --downsample is an integer if --downsample-mode is cell
    if args.downsample and args.downsample_mode in ('cell', 'group'):
        if int(args.downsample) != args.downsample:
            parser.error('`--downsample` must be an integer when using `--downsample-mode cell/group`')
    if args.downsample_mode == 'group' and not args.groups:
        parser.error('`--groups` must be provided when using `--downsample-mode group`')

    from .estimate import estimate
    estimate(
        args.count_dirs,
        args.o,
        args.reads,
        groups=groups,
        ignore_groups_for_pi=args.ignore_groups_for_pi,
        genes=genes,
        downsample=args.downsample,
        downsample_mode=args.downsample_mode,
        cell_threshold=args.cell_threshold,
        cell_gene_threshold=args.cell_gene_threshold,
        control_p_e=control_p_e,
        control=args.control,
        n_threads=args.t,
        temp_dir=temp_dir,
        nasc=args.nasc,
        by_name=args.gene_names,
        seed=args.seed
    )


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'align': parse_align,
    'consensus': parse_consensus,
    'count': parse_count,
    'estimate': parse_estimate,
}


@logger.namespaced('main')
def main():
    parser = argparse.ArgumentParser(
        description=f'{__version__} Complete splicing and labeling quantification from metabolic labeling scRNA-seq'
    )
    parser._actions[0].help = parser._actions[0].help.capitalize()
    parser.add_argument('--list', help='Display list of supported single-cell technologies', action='store_true')
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='<CMD>',
    )

    # Add common options to this parent parser
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('--tmp', metavar='TMP', help='Override default temporary directory', type=str, default='tmp')
    parent.add_argument('--keep-tmp', help='Do not delete the tmp directory', action='store_true')
    parent.add_argument('--verbose', help='Print debugging information', action='store_true')
    parent.add_argument('-t', metavar='THREADS', help='Number of threads to use (default: 8)', type=int, default=8)

    # Command parsers
    parser_ref = setup_ref_args(subparsers, parent)
    parser_align = setup_align_args(subparsers, parent)
    parser_consensus = setup_consensus_args(subparsers, parent)
    parser_count = setup_count_args(subparsers, parent)
    parser_estimate = setup_estimate_args(subparsers, parent)
    command_to_parser = {
        'ref': parser_ref,
        'align': parser_align,
        'consensus': parser_consensus,
        'count': parser_count,
        'estimate': parser_estimate,
    }
    if '--list' in sys.argv:
        print_technologies()

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    logger.debug('Printing verbose output')
    logger.debug(f'Input args: {args}')
    logger.debug(f'Creating {args.tmp} directory')
    if os.path.exists(args.tmp):
        parser.error(
            f'Temporary directory {args.tmp} already exists. '
            'Is another process running? Please specify a different temporary '
            'directory with the `--tmp` option, or remove the one that already '
            'exists.'
        )
    os.makedirs(args.tmp)
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.t)
    try:
        patch_mp_connection_bpo_17560()  # Monkeypatch for python 3.7
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=args.tmp)
        logger.info('Done')
    except Exception:
        logger.exception('An exception occurred')
    finally:
        if not args.keep_tmp:
            logger.debug(f'Removing {args.tmp} directory')
            shutil.rmtree(args.tmp, ignore_errors=True)
