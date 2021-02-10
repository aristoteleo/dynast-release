import argparse
import logging
import os
import shutil
import sys
import warnings

from . import __version__
from .config import RE_CHOICES
from .technology import BARCODE_UMI_TECHNOLOGIES, TECHNOLOGIES_MAP
from .utils import flatten_list

logger = logging.getLogger(__name__)


def silence_logger(name):
    package_logger = logging.getLogger(name)
    package_logger.setLevel(logging.CRITICAL + 10)
    package_logger.propagate = False


def print_technologies():
    """Displays a list of supported technologies along with whether a whitelist
    is provided for that technology.
    """
    headers = ['name', 'whitelist provided', 'barcode (start, length)', 'UMI (start, length)']
    rows = [headers]

    print('List of supported single-cell technologies.\nPositions are 1-indexed.\n')
    for key in sorted(TECHNOLOGIES_MAP):
        t = TECHNOLOGIES_MAP[key]
        barcode = umi = 'NA'
        if t in BARCODE_UMI_TECHNOLOGIES:
            barcode = (t.arguments['--soloCBstart'], t.arguments['--soloCBlen'])
            umi = (t.arguments['--soloUMIstart'], t.arguments['--soloUMIlen'])
        row = [key, 'yes' if t.whitelist_path else '', str(barcode), str(umi)]
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


def setup_ref_args(parser, parent):
    """Helper function to set up a subparser for the `ref` command.

    :param parser: argparse parser to add the `ref` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser

    :return: the newly added parser
    :rtype: argparse.ArgumentParser
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


def setup_align_args(parser, parent):
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
        help='Read strandedness. (default: `forward`)',
        choices=['forward', 'reverse', 'unstranded'],
        default='forward',
    )
    parser_align.add_argument(
        '-w',
        metavar='WHITELIST',
        help=('Path to file of whitelisted barcodes to correct to. '
              'If not provided, all barcodes are used.'),
    )
    parser_align.add_argument('--overwrite', help='Overwrite existing alignment files', action='store_true')
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


def setup_count_args(parser, parent):
    """Helper function to set up a subparser for the `count` command.

    :param parser: argparse parser to add the `count` command to
    :type args: argparse.ArgumentParser
    :param parent: argparse parser parent of the newly added subcommand.
                   used to inherit shared commands/flags
    :type args: argparse.ArgumentParser

    :return: the newly added parser
    :rtype: argparse.ArgumentParser
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
        required=True
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
        help=('BAM tag to use as gene assignments (default: GX)'),
        type=str,
        default='GX',
    )
    parser_count.add_argument(
        '--strand',
        help='Read strandedness. (default: `forward`)',
        choices=['forward', 'reverse', 'unstranded'],
        default='forward',
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
        '--re',
        metavar='RE',
        help=(f'Re-do a step in the pipeline. Available choices are: {", ".join(RE_CHOICES)}.'),
        type=str,
        choices=RE_CHOICES,
        default=None
    )
    parser_count.add_argument(
        '--p-group-by',
        metavar='GROUPBY',
        help=argparse.SUPPRESS,
        type=str,
        default='barcode',
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
        default=False,
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
        help=('Textfile containing filtered cell barcodes. Only these barcodes will '
              'be processed.'),
        type=str,
    )
    parser_count.add_argument(
        '--read-threshold',
        metavar='THRESHOLD',
        help=('Do not attempt statistical correction if there are less than this '
              'many reads. (default: 16)'),
        type=int,
        default=16,
    )
    parser_count.add_argument(
        '--no-velocity',
        '--transcriptome-only',
        help=(
            'Do not prepare matrices for RNA velocity estimation and ignore reads that '
            'are not assigned to the transcriptome.'
        ),
        action='store_true'
    )
    parser_count.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )
    parser_count.add_argument('--seed', help=argparse.SUPPRESS, type=int, default=None)
    parser_count.add_argument(
        '--control',
        help=(
            'Indicate this is a control sample, which is used to estimate the background mutation rate '
            'and/or detect SNPs. The estimated background mutation rate and/or detected SNPs can be '
            'used when running subsequent test samples.'
        ),
        action='store_true',
    )
    parser_count.add_argument(
        '--correct',
        help=(
            'Perform statistical correction of unlabeled and labeled read counts. '
            'This option can be used multiple times to correct multiple species. '
            'By default, no correction is performed.'
        ),
        action='append',
        default=None,
        choices=['total', 'transcriptome', 'spliced', 'unspliced']
    )
    parser_count.add_argument(
        '--p-e',
        help='Textfile containing a single number, indicating the estimated background mutation rate',
        type=str,
        default=None,
    )
    parser_count.add_argument('--subset-threshold', help=argparse.SUPPRESS, type=int, default=8000)
    parser_count.add_argument(
        'bam',
        help=(
            'Alignment BAM file that contains the appropriate UMI and barcode tags, '
            'specifiable with `--umi-tag`, and `--barcode-tag`.'
        ),
        type=str,
    )

    return parser_count


def parse_ref(parser, args, temp_dir=None):
    """Parser for the `ref` command.
    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
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

    from .align import align
    align(
        args.fastqs,
        args.i,
        args.o,
        technology,
        whitelist_path=args.w,
        n_threads=args.t,
        temp_dir=temp_dir,
        nasc=args.nasc,
        overrides=overrides
    )


def parse_count(parser, args, temp_dir=None):
    """Parser for the `count` command.
    :param args: Command-line arguments dictionary, as parsed by argparse
    :type args: dict
    """
    # Check quality
    if args.quality < 0 or args.quality > 41:
        parser.error('`--quality` must be in [0, 42)')

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

    # Check group by
    if args.p_group_by.lower() == 'none':
        args.p_group_by = None
    else:
        args.p_group_by = args.p_group_by.split(',')

    # Read barcodes
    barcodes = None
    if args.barcodes:
        with open(args.barcodes, 'r') as f:
            barcodes = [line.strip() for line in f if not line.isspace()]
        logger.info(f'Ignoring cell barcodes not in the {len(barcodes)} barcodes provided by `--barcodes`')
    else:
        logger.warning('`--barcodes` not provided. All cell barcodes will be processed.')

    if not args.correct:
        logger.warning('No statistical correction will be performed because `--correct` is not provided.')
    elif 'total' in args.correct and args.no_velocity:
        parser.error(
            '`--no-velocity` or `--transcriptome-only` can not be used with `--correct total`. '
            'Use `--correct transcriptome` instead.'
        )

    if not args.gene_tag:
        parser.error('`--gene-tag` must be a valid tag')

    # Check conversions
    if args.nasc and args.conversion != ['TC']:
        parser.error('`--conversion TC` is required for `--nasc`')
    conversions = []
    for arg in args.conversion:
        conversions.append(arg.upper().split(','))
    flattened = list(flatten_list(conversions))
    if len(set(flattened)) != len(flattened):
        parser.error('duplicate conversions are not allowed for `--conversion`')
    if set(f[0] for f in flattened) == {'A', 'C', 'G', 'T'} and control_p_e is None:
        parser.error('`--p-e` must be specified when conversions apply to all four nucleotides')

    from .count import count
    count(
        args.bam,
        args.g,
        args.o,
        strand=args.strand,
        umi_tag=args.umi_tag,
        barcode_tag=args.barcode_tag,
        gene_tag=args.gene_tag,
        barcodes=barcodes,
        control=args.control,
        quality=args.quality,
        conversions=conversions,
        snp_threshold=args.snp_threshold,
        snp_csv=args.snp_csv,
        correct=args.correct,
        read_threshold=args.read_threshold,
        control_p_e=control_p_e,
        p_group_by=args.p_group_by,
        n_threads=args.t,
        re=args.re,
        temp_dir=temp_dir,
        nasc=args.nasc,
        subset_threshold=args.subset_threshold,
        velocity=not args.no_velocity,
        seed=args.seed,
    )


COMMAND_TO_FUNCTION = {
    'ref': parse_ref,
    'align': parse_align,
    'count': parse_count,
}


def main():
    parser = argparse.ArgumentParser(description=f'{__version__}')
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
    parser_count = setup_count_args(subparsers, parent)
    command_to_parser = {
        'ref': parser_ref,
        'align': parser_align,
        'count': parser_count,
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

    logging.basicConfig(
        format='[%(asctime)s] %(levelname)7s %(message)s',
        level=logging.DEBUG if args.verbose else logging.INFO,
        force=True,
    )

    # Silence logging from other packages
    silence_logger('anndata')
    silence_logger('h5py')
    silence_logger('numba')
    silence_logger('pysam')
    silence_logger('pystan')

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
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=args.tmp)
    except Exception:
        logger.exception('An exception occurred')
    finally:
        if not args.keep_tmp:
            logger.debug(f'Removing {args.tmp} directory')
            shutil.rmtree(args.tmp, ignore_errors=True)
